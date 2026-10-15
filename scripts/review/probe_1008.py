"""Inspect schemas and survey vocabulary; no participant-level output."""

from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.functions.paths import PROJECT_ROOT

import json
import os
from google.cloud import bigquery

ROOT = PROJECT_ROOT


def main():
    client = bigquery.Client(project=os.environ['GOOGLE_PROJECT'])
    cdr = os.environ['WORKSPACE_CDR']
    out = {}
    for table in [
        'person',
        'cb_search_person',
        'ds_survey',
        'observation',
        'death',
        'aou_death',
        'procedure_occurrence',
    ]:
        out[table] = [
            {"name": f.name, "type": f.field_type}
            for f in client.get_table(f'{cdr}.{table}').schema
        ]
    sql = f"""SELECT concept_id, concept_name FROM `{cdr}.concept`
    WHERE REGEXP_CONTAINS(LOWER(concept_name), r'family.*(breast|cancer)|breast.*(family|mother|father|sibling)|personal and family|family health history|mastectomy')
    ORDER BY concept_id"""
    cfg = bigquery.QueryJobConfig(maximum_bytes_billed=40 * 1024**3)
    out['concepts'] = client.query(sql, job_config=cfg).to_dataframe().to_dict('records')
    sql = f"""SELECT DISTINCT survey, survey_version_concept_id, survey_version_name,
    question_concept_id, question, answer_concept_id, answer
    FROM `{cdr}.ds_survey`
    WHERE REGEXP_CONTAINS(LOWER(survey), r'family|personal.*health')
    AND REGEXP_CONTAINS(LOWER(question), r'cancer|relative|family')"""
    out['survey_templates'] = client.query(sql, job_config=cfg).to_dataframe().to_dict('records')
    sql = f"""SELECT concept_id, concept_name, concept_code, vocabulary_id, concept_class_id
    FROM `{cdr}.concept` WHERE vocabulary_id IN ('PPI', 'AoU')
    AND REGEXP_CONTAINS(LOWER(CONCAT(concept_name,' ',concept_code)), r'breast|cancer|consent|enroll')
    ORDER BY concept_id"""
    out['source_concepts'] = client.query(sql, job_config=cfg).to_dataframe().to_dict('records')
    sql = f"""SELECT DISTINCT o.observation_concept_id,o.observation_source_concept_id,
    c.concept_name,c.concept_code,c.vocabulary_id
    FROM `{cdr}.observation` o
    LEFT JOIN `{cdr}.concept` c ON o.observation_source_concept_id=c.concept_id
    WHERE o.observation_concept_id IN (836772,43528515)
    ORDER BY o.observation_concept_id,o.observation_source_concept_id"""
    cfg = bigquery.QueryJobConfig(maximum_bytes_billed=100 * 1024**3)
    out['observed_survey_source_templates'] = (
        client.query(sql, job_config=cfg).to_dataframe().to_dict('records')
    )
    path = ROOT / 'results/report/62_review_1008_schema_probe.json'
    path.write_text(json.dumps(out, indent=2, default=int) + '\n')
    print('Schema and public vocabulary probe saved; no participant rows printed.')


if __name__ == '__main__':
    main()
