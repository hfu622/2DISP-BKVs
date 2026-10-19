"""Adversarial checks: a solver's success flag cannot certify invalid geometry."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from verify_layout import verify

ITEM = {'id':0,'demand':2,'orientation':{'rotation':{'mode':'discrete','angles':[0]}},
        'shape':{'type':'simple_polygon','data':[[0,0],[1,0],[1,1],[0,1]]}}
INPUT = {'name':'fixture','strip_height':1,'items':[ITEM]}
MANIFEST = {'instance':'fixture','fixed_width':1,'total_area':2,'xml_sha256':'fixture', 'boundary_discrepancies':[]}
OUTPUT = {'instance':INPUT, 'solution':{'strip_width':2, 'density':1,
    'layout':{'placed_items':[{'item_id':0,'transformation':{'rotation':0,'translation':[0,0]}},
                             {'item_id':0,'transformation':{'rotation':0,'translation':[1,0]}}]}}}


class VerificationTests(unittest.TestCase):
    def check_output(self, raw, **kwargs):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'solution.json'
            path.write_text(json.dumps(raw),encoding='utf-8')
            with patch('verify_layout.load_xml', return_value=(INPUT,MANIFEST)):
                return verify('unused.xml',path, **kwargs)

    def test_declared_tiny_width_is_allowed_and_both_densities_reported(self):
        raw=copy.deepcopy(OUTPUT);raw['instance']['strip_height']=1.0001
        check=self.check_output(raw, allowed_width_expansion_ratio=0.0001)
        self.assertTrue(check['passed_with_declared_tolerances'])
        self.assertAlmostEqual(check['utilization_percent_actual_width'],100/1.0001)
        self.assertEqual(check['utilization_percent_recomputed'],100)

    def test_large_expansion_and_shrinking_rejected(self):
        for width in [1.001,0.9999]:
            raw=copy.deepcopy(OUTPUT);raw['instance']['strip_height']=width
            check=self.check_output(raw, allowed_width_expansion_ratio=0.0001)
            self.assertFalse(check['passed_with_declared_tolerances'])

    def test_exact_contact_allowed_and_density_recomputed(self):
        raw=copy.deepcopy(OUTPUT);raw['solution']['density']=0.123
        check=self.check_output(raw)
        self.assertTrue(check['passed_with_declared_tolerances'])
        self.assertTrue(check['zero_area_violation'])
        self.assertEqual(check['utilization_percent_recomputed'],100)

    def test_overlap_rejected(self):
        raw=copy.deepcopy(OUTPUT);raw['solution']['layout']['placed_items'][1]['transformation']['translation']=[0.5,0]
        self.assertIn('positive-area overlap exceeds tolerance',self.check_output(raw)['errors'])

    def test_boundary_violation_rejected(self):
        raw=copy.deepcopy(OUTPUT);raw['solution']['layout']['placed_items'][1]['transformation']['translation']=[1.01,0]
        self.assertFalse(self.check_output(raw)['passed_with_declared_tolerances'])

    def test_missing_copy_rejected(self):
        raw=copy.deepcopy(OUTPUT);raw['solution']['layout']['placed_items'].pop()
        self.assertIn('missing / duplicate / unknown items',self.check_output(raw)['errors'])

    def test_forbidden_rotation_rejected(self):
        raw=copy.deepcopy(OUTPUT);raw['solution']['layout']['placed_items'][0]['transformation']['rotation']=0.1
        self.assertIn('forbidden rotation',self.check_output(raw)['errors'])

    def test_full_turn_representatives_are_equivalent_but_added_angles_rejected(self):
        raw=copy.deepcopy(OUTPUT);raw['instance']['items'][0]['orientation']['rotation']['angles']=[-360]
        self.assertTrue(self.check_output(raw)['passed_with_declared_tolerances'])
        raw['instance']['items'][0]['orientation']['rotation']['angles']=[-360,90]
        self.assertIn('input allowed rotations changed',self.check_output(raw)['errors'])

    def test_nonfinite_transform_rejected(self):
        raw=copy.deepcopy(OUTPUT);raw['solution']['layout']['placed_items'][0]['transformation']['translation'][0]=float('nan')
        self.assertIn('nonfinite transformation',self.check_output(raw)['errors'])

    def test_casefold_name_allowed_but_different_instance_rejected(self):
        raw=copy.deepcopy(OUTPUT);raw['instance']['name']='FIXTURE'
        self.assertTrue(self.check_output(raw)['passed_with_declared_tolerances'])
        raw['instance']['name']='other'
        self.assertIn('instance name changed',self.check_output(raw)['errors'])


if __name__=='__main__':
    unittest.main()
