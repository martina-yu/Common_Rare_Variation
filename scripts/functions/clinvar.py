# Documentation: scripts/README.md.
import hashlib
import gzip
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import time
from urllib.parse import unquote
import xml.etree.ElementTree as ET

import pandas as pd
import requests
import yaml

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'data/29_rowlands_extension/section6'
XML = OUT / 'clinvar_xml'
PLP = {'pathogenic', 'likely pathogenic', 'pathogenic/likely pathogenic'}
STARS = {
    'criteria provided, multiple submitters, no conflicts': 2,
    'criteria provided, single submitter': 1,
    'reviewed by expert panel': 3,
    'practice guideline': 4,
}


def write_tsv_if_changed(frame, path):
    """Keep valid genotype checkpoints when a reference revalidation is identical."""
    text = frame.to_csv(index=False, sep='\t')
    if path.exists():
        with gzip.open(path, 'rt') as f:
            if f.read() == text:
                return
    frame.to_csv(path, sep='\t', index=False, compression={'method': 'gzip', 'mtime': 0})


def clean(text):
    return re.sub(r'\s+', ' ', unquote(text or '').replace('_', ' ')).strip().lower()


def exclusions(values):
    text = ' | '.join(clean(v) for v in values)
    if 'low penetrance' in text:
        return 'low_penetrance'
    if re.search(r'\brisk (?:factor|allele)\b', text):
        return 'risk_factor_or_risk_allele'
    conflict_text = re.sub(r'\bno conflicts?\b', '', text)
    if 'conflict' in conflict_text or 'benign' in text:
        return 'conflict_or_opposing_B_LB'
    return ''


def fetch(ids):
    for attempt in range(6):
        try:
            params = {'db': 'clinvar', 'id': ','.join(ids), 'rettype': 'vcv'}
            if all(x.isdigit() for x in ids):
                params['is_variationid'] = 'true'
            r = requests.get(
                'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',
                params=params,
                timeout=90,
            )
            r.raise_for_status()
            root = ET.fromstring(r.content)
            records = root.findall('.//VariationArchive')
            if not records:
                raise ValueError('No VariationArchive returned')
            time.sleep(0.4)
            return records
        except (requests.RequestException, ET.ParseError, ValueError):
            if attempt == 5:
                raise RuntimeError('Public ClinVar XML retrieval failed') from None
            time.sleep(min(2**attempt, 15))


def historical_xml(ids, cutoff, hold_unavailable=False):
    XML.mkdir(parents=True, exist_ok=True)
    missing = [v for v in ids if not (XML / f'{v}.xml').exists()]
    for start in range(0, len(missing), 20):
        batch = missing[start : start + 20]
        records = {v.attrib['VariationID']: v for v in fetch(batch)}
        assert set(records) == set(batch), 'ClinVar did not return every requested ID'
        for vid, v in records.items():
            version = int(v.attrib['Version'])
            while v.attrib['DateLastUpdated'] > cutoff:
                version -= 1
                if version <= 0 and hold_unavailable:
                    v = None
                    break
                assert version > 0, 'No historical VCV before the frozen release'
                accession = f"{v.attrib['Accession']}.{version}"
                v = fetch([accession])[0]
                assert (
                    int(v.attrib['Version']) == version
                ), 'Historical endpoint ignored the version'
            if v is not None:
                ET.ElementTree(v).write(XML / f'{vid}.xml', encoding='utf-8', xml_declaration=True)
        print(f'Historical ClinVar XML batch {start // 20 + 1} completed.', flush=True)
    result = {}
    for vid in ids:
        if hold_unavailable and not (XML / f'{vid}.xml').exists():
            continue
        v = ET.parse(XML / f'{vid}.xml').getroot()
        assert v.attrib['VariationID'] == vid and v.attrib['DateLastUpdated'] <= cutoff
        result[vid] = v
    return result


def disease_rule(gene, name, inheritance=''):
    """Disease rule."""
    name, inheritance = clean(name), clean(inheritance)
    if 'recessive' in inheritance or 'biallelic' in inheritance:
        return False, 'assertion_supports_recessive_or_biallelic_disease'
    if any(
        x in name
        for x in [
            'fanconi',
            'ataxia-telangiectasia',
            'ataxia telangiectasia',
            'seckel',
            'developmental delay',
            'autism',
            'cleft lip',
            'blepharocheilodontic',
        ]
    ):
        return False, 'unrelated_or_biallelic_disease'
    common = {
        f'{gene.lower()}-related cancer predisposition',
        f'{gene.lower()}-related cancer susceptibility',
    }
    accepted = name in common
    if gene in {'BRCA1', 'BRCA2'}:
        accepted |= any(
            x in name
            for x in [
                'hereditary breast and ovarian cancer',
                'hereditary breast ovarian cancer',
                'breast-ovarian cancer, familial',
                'breast and/or ovarian cancer',
            ]
        )
    elif gene == 'TP53':
        accepted |= 'li-fraumeni' in name or 'li fraumeni' in name
    elif gene == 'PTEN':
        accepted |= any(
            x in name
            for x in [
                'pten hamartoma',
                'cowden',
                'bannayan-riley-ruvalcaba',
                'bannayan riley ruvalcaba',
            ]
        )
    elif gene == 'CDH1':
        accepted |= any(
            x in name
            for x in [
                'hereditary diffuse gastric',
                'diffuse gastric and lobular breast',
                'gastric cancer, hereditary diffuse',
                'hereditary lobular breast',
            ]
        )
    elif gene == 'NF1':
        accepted |= name in {
            'neurofibromatosis, type 1',
            'neurofibromatosis type 1',
            'neurofibromatosis 1',
            'von recklinghausen disease',
        }
    elif gene == 'STK11':
        accepted |= 'peutz-jeghers' in name or 'peutz jeghers' in name
    # Named breast cancer susceptibility assertions are relevant only for the
    # prespecified breast-cancer panel and absent a recessive assertion.
    if gene in {'BRCA1', 'BRCA2', 'PALB2', 'CHEK2', 'ATM', 'BARD1', 'RAD51C', 'RAD51D'}:
        accepted |= name in {
            'breast cancer',
            'hereditary breast cancer',
            'familial breast cancer',
            'breast cancer, susceptibility to',
            'hereditary breast cancer susceptibility',
        }
        accepted |= bool(re.fullmatch(r'breast cancer, (?:familial|susceptibility to),? \d+', name))
        accepted |= name.startswith('hereditary breast and ovarian cancer')
    if accepted:
        if inheritance and 'dominant' not in inheritance:
            return False, 'inheritance_attribute_requires_review'
        return True, 'AD_named_syndrome_or_monoallelic_breast_cancer_susceptibility'
    return False, 'unmapped_disease_requires_review'


def review_variant(v, gene, frozen_sig, frozen_review, min_stars=2):
    if min_stars not in (1, 2):
        raise ValueError('Supported prespecified ClinVar thresholds are 1 and 2')
    rows, eligible_rcvs = [], []
    contradictory_inheritance = False
    # Include all current germline submissions, even noncontributing ones, so
    # an expert-panel aggregate does not silently hide opposing B/LB evidence.
    assertions = []
    for scv in v.findall('.//ClinicalAssertion'):
        if scv.findtext('RecordStatus', 'current') != 'current':
            continue
        sig = scv.findtext('Classification/GermlineClassification', '')
        if not sig:
            continue
        acc = scv.find('ClinVarAccession')
        diseases = [
            n.text or ''
            for n in scv.findall('TraitSet/Trait/Name/ElementValue')
            if n.attrib.get('Type') == 'Preferred'
        ]
        inheritance = ';'.join(
            n.text or ''
            for n in scv.findall('.//Attribute')
            if 'inheritance' in n.attrib.get('Type', '').lower()
        )
        assertions.append(sig)
        for disease in diseases or ['']:
            mode = clean(inheritance)
            if (
                clean(sig) in PLP
                and disease_rule(gene, disease, '')[0]
                and ('recessive' in mode or 'biallelic' in mode)
                and 'dominant' not in mode
            ):
                contradictory_inheritance = True
            ok, reason = disease_rule(gene, disease, inheritance)
            rows.append(
                dict(
                    record_type='SCV',
                    accession=acc.attrib['Accession'] if acc is not None else '',
                    version=acc.attrib.get('Version', '') if acc is not None else '',
                    classification=sig,
                    review_status=scv.findtext('Classification/ReviewStatus', ''),
                    disease=disease,
                    disease_ids='|'.join(
                        f"{x.attrib.get('DB', '')}:{x.attrib.get('ID', '')}"
                        for x in scv.findall('TraitSet/Trait/XRef')
                    ),
                    inheritance=inheritance,
                    disease_relevant=ok,
                    mapping_reason=reason,
                    contributes=scv.attrib.get('ContributesToAggregateClassification', ''),
                    germline_origin=';'.join(x.text or '' for x in scv.findall('.//Sample/Origin')),
                )
            )
    block = exclusions(assertions)
    aggregate = v.find('ClassifiedRecord/Classifications/GermlineClassification')
    if aggregate is not None:
        xml_sig = aggregate.findtext('Description', '')
        xml_review = aggregate.findtext('ReviewStatus', '')
    else:
        xml_sig, xml_review = '', ''
    snapshot_match = clean(xml_sig) == clean(frozen_sig) and clean(xml_review) == clean(
        frozen_review
    )
    for rcv in v.findall('.//RCVAccession'):
        cl = rcv.find('RCVClassifications/GermlineClassification')
        if cl is None:
            continue
        sig, rev = cl.findtext('Description', ''), cl.findtext('ReviewStatus', '')
        names = [x.text or '' for x in rcv.findall('ClassifiedConditionList/ClassifiedCondition')]
        mapped = [disease_rule(gene, x) for x in names]
        relevant = bool(mapped) and all(ok for ok, _ in mapped)
        rows.append(
            dict(
                record_type='RCV',
                accession=rcv.attrib['Accession'],
                version=rcv.attrib.get('Version', ''),
                classification=sig,
                review_status=rev,
                disease='|'.join(names),
                disease_ids='|'.join(
                    f"{x.attrib.get('DB', '')}:{x.attrib.get('ID', '')}"
                    for x in rcv.findall('ClassifiedConditionList/ClassifiedCondition')
                ),
                inheritance='AD (curated syndrome mapping)' if relevant else '',
                disease_relevant=relevant,
                mapping_reason=';'.join(reason for _, reason in mapped),
                contributes='',
                germline_origin='germline classification field',
            )
        )
        if relevant and clean(sig) in PLP and STARS.get(clean(rev), 0) >= min_stars:
            eligible_rcvs.append(rcv.attrib['Accession'])
    # A single unrelated benign assertion triggers the user-requested blanket
    # exclusion, while disease suitability still depends on relevant P/LP RCVs.
    block = block or exclusions([r['classification'] for r in rows])
    status = block or ('historical_XML_VCF_aggregate_mismatch' if not snapshot_match else '')
    if not status and contradictory_inheritance:
        status = 'explicit_recessive_inheritance_contradicts_AD_disease_mapping_pending_review'
    if not status and not eligible_rcvs:
        status = f'no_relevant_at_least_{min_stars}star_P_LP_RCV_pending_review'
    return not status, status or 'eligible', eligible_rcvs, rows
