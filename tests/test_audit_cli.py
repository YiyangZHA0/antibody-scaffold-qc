"""CLI safety/routing regressions, independent of numbering-model execution."""
import json
from pathlib import Path
import pytest
from scaffold.cli import parser, execute
from scaffold.errors import ScaffoldError


def report(path, **kwargs):
    if path.name=='bad.fasta':
        raise ScaffoldError('invalid test input')
    return {'input':str(path),'status':'not_assessable','domain_count':0,'domains':[],
            'risks':[],'risk_counts':{},'limitations':[],'coordinate_provenance':None,'structure_warnings':[]}


def test_mixed_jobs_continue_and_do_not_overwrite_inputs(tmp_path,monkeypatch):
    from scaffold import audit_cli
    monkeypatch.setattr(audit_cli,'audit_file',lambda path,*a,**kw:report(path))
    files=[tmp_path/name for name in ['good.pdb','bad.fasta','good.cif']]
    for p in files:p.write_text('unchanged input')
    out=tmp_path/'results'
    args=parser().parse_args(['audit',*map(str,files),'--out',str(out)])
    result=execute(args,numberer=object())
    assert result['files']==3 and result['failed_files']==1 and not result['passed']
    assert len({r['output'] for r in result['results']})==3
    assert all(p.read_text()=='unchanged input' for p in files)
    assert (out/'summary.json').is_file()
    assert sum((Path(r['output'])/'audit.json').is_file() for r in result['results'])==2


def test_output_directory_is_never_reused(tmp_path,monkeypatch):
    p=tmp_path/'input.fasta';p.write_text('>a\nAAAA')
    out=tmp_path/'result';out.mkdir();(out/'keep.txt').write_text('keep')
    args=parser().parse_args(['audit',str(p),'--out',str(out)])
    with pytest.raises(ScaffoldError,match='not empty'):
        execute(args,numberer=object())
    assert (out/'keep.txt').read_text()=='keep'


@pytest.mark.parametrize('option',['--sequence','--reference'])
def test_paired_evidence_cannot_be_silently_reused_across_jobs(tmp_path,option):
    files=[tmp_path/name for name in ['a.pdb','b.pdb']]
    for p in files:p.touch()
    args=parser().parse_args(['audit',*map(str,files),option,'reference.fasta'])
    with pytest.raises(ScaffoldError,match='one input only'):
        execute(args,numberer=object())


def test_fasta_chain_options_rejected(tmp_path):
    p=tmp_path/'input.fasta';p.write_text('>a\nAAAA')
    args=parser().parse_args(['audit',str(p),'--chain','A'])
    with pytest.raises(ScaffoldError,match='Chain options'):
        execute(args,numberer=object())
