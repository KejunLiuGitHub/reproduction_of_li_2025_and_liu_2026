#!/usr/bin/env python3
"""Deterministic geometry/graph tests; synthetic arrays are not scientific inputs."""
import argparse
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import unittest
import numpy as np

sys.dont_write_bytecode = True
SOURCE = Path(__file__).with_name('analyze_network.py')
spec = importlib.util.spec_from_file_location('network', SOURCE)
n = importlib.util.module_from_spec(spec)
spec.loader.exec_module(n)
P = n.ROOT / 'campaign/li2025_recovery_20261111/inputs/p_interface.xyz'


def synthetic():
    x = np.full((296, 3), 800.0)
    for w in range(96):
        x[3*w+2] = [100 + 5*w, 200, 200]
        x[3*w] = x[3*w+2] + [1, 0, 0]
        x[3*w+1] = x[3*w+2] + [0, 1, 0]
    x[2], x[0], x[1], x[289] = [10,10,10], [11,10,10], [10,11,10], [13,10,10]
    return x


def identities(x, box, criterion='standard'):
    return {e[:3] for e in n.edges_from_geometry(n.bond_geometry(x, box), n.THRESHOLDS[criterion])}


class NetworkTests(unittest.TestCase):
    def test_known_linear_pair_and_self_exclusion(self):
        edges = n.edges_from_geometry(n.bond_geometry(synthetic(), np.array([1000.]*3)), (2.5,150))
        match = [e for e in edges if e[:3] == (3,1,290)]
        self.assertEqual(len(match), 1)
        self.assertAlmostEqual(match[0][3], 2.0)
        self.assertAlmostEqual(match[0][4], 180)
        self.assertTrue(all(e[0] != e[2] for e in edges))

    def test_distance_and_angle_sensitivity(self):
        x, box = synthetic(), np.array([1000.]*3)
        x[289] = [13.45,10,10]
        self.assertNotIn((3,1,290), identities(x, box, 'strict'))
        self.assertIn((3,1,290), identities(x, box))
        x[289] = [13.55,10,10]
        self.assertNotIn((3,1,290), identities(x, box))
        self.assertIn((3,1,290), identities(x, box, 'loose'))
        for angle, expect in [(149, False), (151, True)]:
            theta = np.radians(180-angle)
            x[289] = x[0] + 2*np.array([np.cos(theta),np.sin(theta),0])
            self.assertEqual((3,1,290) in identities(x,box), expect)

    def test_periodic_boundary_pair(self):
        x = synthetic()
        x[2], x[0], x[1], x[289] = [13,10,10], [14,10,10], [13,11,10], [.2,10,10]
        edges = n.edges_from_geometry(n.bond_geometry(x,np.array([15.,15.,30.])),(2.5,150))
        e = next(e for e in edges if e[:3] == (3,1,290))
        self.assertAlmostEqual(e[3],1.2)
        self.assertAlmostEqual(e[4],180)

    def test_whole_box_individual_translation(self):
        _, x = next(n.sol.frames(P))
        box = np.array([15.,15.,30.])
        shifts = np.random.default_rng(33).integers(-3,4,(296,3))*box
        for name in n.THRESHOLDS:
            self.assertEqual(identities(x,box,name),identities(x+shifts+[.3,-.4,.7],box,name))

    def test_xy_rotation_preserves_orthorhombic_cell(self):
        _, x = next(n.sol.frames(P))
        box = np.array([15.,15.,30.])
        rotation = np.array([[0,-1,0],[1,0,0],[0,0,1]])
        for name in n.THRESHOLDS:
            self.assertEqual(identities(x,box,name),identities(x@rotation.T,box,name))

    def test_registered_counts_and_donor_identity(self):
        source = n.ROOT / 'results/li2025_autofollow/20261113T213847+0800/analysis/p_iface_v2/SOLVATION.json'
        ref = json.loads(source.read_text())
        for i,(step,x) in enumerate(n.sol.frames(ref['source'])):
            if i not in (0,1,154,226,230,231,500):
                continue
            old = n.sol.describe(x,np.array(ref['box']),interface=False)
            geometry = n.bond_geometry(x,np.array(ref['box']))
            self.assertEqual(int(np.sum(geometry[3]>1.3)),old['water_parent_assignment_warning'])
            for name,threshold in n.THRESHOLDS.items():
                edges = n.edges_from_geometry(geometry,threshold)
                graph = n.local_graph(edges)
                suffix = '' if name == 'standard' else '_'+name
                for role,atom in n.ROLES.items():
                    self.assertEqual(graph['solute_HB_counts'][role],old['hb_'+role+suffix])
                    if name == 'standard':
                        self.assertEqual(sorted(e[1] for e in edges if e[2]==atom),old['donor_H_1based'][role])

    def test_parent_oh_warning_is_not_reassignment(self):
        x = synthetic()
        x[0] = x[2] + [1.4,0,0]
        geometry = n.bond_geometry(x,np.array([1000.]*3))
        flagged = (n.sol.WATER_H[geometry[3]>1.3]+1).tolist()
        self.assertEqual(flagged,[1])
        self.assertEqual(n.OWNERS[0],2)

    def test_local_shell_direction_and_family(self):
        edges = [(3,1,290,2,180),(6,4,294,2,180),(3,2,9,2,180),(9,7,6,2,180),(9,8,12,2,180)]
        g = n.local_graph(edges)
        self.assertEqual(g['first_shell'],[3,6])
        self.assertEqual(g['second_shell'],[9])
        self.assertEqual(g['min_water_path_undirected'],2)
        self.assertEqual(g['min_water_path_Oc_donor_waters_to_Os_donor_waters'],2)
        self.assertIsNone(g['min_water_path_Os_donor_waters_to_Oc_donor_waters'])
        self.assertEqual(g['water_degree_families_global']['9'],'D2A1')
        self.assertNotIn([9,8,12],g['edges_1based_D_H_A'])

    def test_continuous_intermittent_and_censoring(self):
        steps=list(range(0,70,10));idx=[0,1,3,6]
        a=n.episodes(idx,steps);b=n.episodes(idx,steps,1)
        self.assertEqual(len(a),3);self.assertEqual(len(b),2)
        self.assertEqual(a[0]['endpoint_span_fs'],10)
        self.assertEqual(b[0]['endpoint_span_fs'],30)
        self.assertEqual(b[0]['absent_frames_within_span'],1)
        self.assertTrue(a[0]['left_censored']);self.assertTrue(a[-1]['right_censored'])
        self.assertEqual(a[-1]['endpoint_span_fs'],0)


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path);a=parser.parse_args()
    stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(NetworkTests))
    print(stream.getvalue())
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps({'status':'PASS' if result.wasSuccessful() else 'FAIL',
            'tests':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),'details':stream.getvalue(),
            'source_sha256':{str(p.resolve()):n.sha(p) for p in [SOURCE,Path(__file__),n.SOL,P]}},indent=2)+'\n')
    sys.exit(0 if result.wasSuccessful() else 1)
