"""Original ALT-index-aware QC and in-memory normalization."""
import math
import subprocess
from pathlib import Path

def call_qc(gt, ad, gq, target_index, n_alts, gene):
    try:
        alleles = [int(x) for x in gt.replace('|', '/').split('/')]
        if len(alleles) != 2 or any((x < 0 or x > n_alts for x in alleles)):
            return None
        dose = alleles.count(target_index)
        if not dose:
            return None
        gq = float(gq)
        depths = [int(x) for x in ad.split(',')]
        if len(depths) != n_alts + 1 or any((x < 0 for x in depths)) or sum(depths) <= 0:
            return {'gq_pass': math.isfinite(gq) and gq >= 20, 'valid_AD': False}
        vaf = depths[target_index] / sum(depths)
        main = vaf >= 0.8 if dose == 2 else 0.3 <= vaf <= 0.8
        sensitive = vaf >= 0.8 if dose == 2 else 0.3 <= vaf <= (0.7 if gene in {'TP53', 'NF1'} else 0.8)
        return dict(gq_pass=math.isfinite(gq) and gq >= 20, valid_AD=True, AD_total=sum(depths), vaf=vaf, dose=dose, main=main, sensitive=sensitive)
    except (ValueError, IndexError):
        try:
            alleles = [int(x) for x in gt.replace('|', '/').split('/')]
            if len(alleles) == 2 and target_index in alleles:
                return {'gq_pass': False, 'valid_AD': False}
        except ValueError:
            pass
        return None

def normalized_native_sites(raw, fasta):
    """Normalize single ALT site copies; leave original GT/AD indices untouched."""
    tokens, lines = ({}, ['##fileformat=VCFv4.2'])
    for l in Path(str(fasta) + '.fai').read_text().splitlines():
        chrom, length, *_ = l.split('\t')
        lines.append(f'##contig=<ID={chrom},length={length}>')
    lines.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
    for native in raw:
        chrom, pos, ref, alts, filt = native
        for index, alt in enumerate(alts.split(','), 1):
            if any((x not in 'ACGTN' for x in alt)):
                continue
            token = f'n{len(tokens)}'
            tokens[token] = (chrom, pos, ref, alts, index, filt)
            lines.append(f'{chrom}\t{pos}\t{token}\t{ref}\t{alt}\t.\t.\t.')
    if not tokens:
        return {}
    norm = subprocess.run(['bcftools', 'norm', '-f', str(fasta), '-c', 'e', '-m', '-any', '-Ov', '-'], input='\n'.join(lines) + '\n', capture_output=True, text=True, check=True)
    mapping = {}
    for l in norm.stdout.splitlines():
        if l.startswith('#'):
            continue
        chrom, pos, token, ref, alt, *_ = l.split('\t')
        mapping[tokens[token]] = f'{chrom}:{pos}:{ref}:{alt}'
    assert len(mapping) == len(tokens)
    return mapping


def source_site_pass(filters):
    """Missing site information fails; an empty defined set is PASS."""
    return filters is not None and (len(filters) == 0 or set(filters) == {'PASS'})


def passing_call(gt, ad, gq, target_index, n_alts, gene, ft, site_filters):
    """DP is not an additional main threshold; AD denominator includes every ALT."""
    if target_index < 1 or target_index > n_alts:
        raise ValueError('Invalid original ALT index')
    if not isinstance(gt, str) or not isinstance(ad, str):
        return {'main': False, 'strict_vaf': False}
    qc = call_qc(gt, ad, gq, target_index, n_alts, gene)
    base = source_site_pass(site_filters) and ft in {'PASS', '.'}
    valid = bool(qc and qc.get('gq_pass') and qc.get('valid_AD'))
    return {'main': bool(base and valid and qc.get('main')),
            'strict_vaf': bool(base and valid and qc.get('sensitive'))}


def local_to_global_gt(local_gt, local_alleles):
    """Convert VDS LGT through LA before inspecting original global ALT dosage."""
    parts = local_gt.replace('|', '/').split('/')
    if len(parts) != 2 or any(p == '.' for p in parts):
        return './.'
    try:
        if any(int(p) < 0 for p in parts):
            raise ValueError
        return '/'.join(str(local_alleles[int(p)]) for p in parts)
    except (ValueError, IndexError):
        raise ValueError('Invalid local genotype mapping') from None
