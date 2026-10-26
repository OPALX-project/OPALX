"""Shared snapshot-extractor block for OPALX Catalyst state scripts.

Exec'd by catalyst-muon-{test-1,fields,combined}.py with a namespace
containing `source` (the producer to snapshot), `name` (its channel
name), and `CreateExtractor` (from paraview.simple).

OPALX_CATALYST_DUMP_SNAPSHOT=N creates a VTPC data extractor that
writes on iteration N and every Nth iteration after that. Set to 0
or unset to disable.
"""
import os

_snapshot_step = 0
try:
    _snapshot_step = int(os.environ.get('OPALX_CATALYST_DUMP_SNAPSHOT', '0'))
except ValueError:
    _snapshot_step = 0

if _snapshot_step > 0:
    ext = CreateExtractor('VTPC', source, registrationName='Snapshot_' + name)
    ext.Trigger = 'Time Step'
    ext.Trigger.Frequency = _snapshot_step
    ext.Writer.FileName = name + '_snapshot/{timestep:06d}.vtpc'
