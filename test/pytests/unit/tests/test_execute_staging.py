"""Tests for ISSM 2026.3-compatible execution staging."""

from pathlib import Path
from types import SimpleNamespace
import tarfile

import pytest

from pyissm.model import execute
from pyissm.model.classes.cluster import _create_input_archive
from pyissm.model.classes.issmsettings import issmsettings
from pyissm.tools import config


class _FakeToolkits:
    def write_toolkits_file(self, filename):
        Path(filename).write_text('toolkits', encoding = 'utf-8')


class _FakeCluster:
    def __init__(self, name, executionpath):
        self.name = name
        self.executionpath = str(executionpath)
        self.uploads = []
        self.launches = []
        self.queue_options = None

    def build_queue_script(self, **kwargs):
        self.queue_options = kwargs
        Path(kwargs['file_prefix'] + '.queue').write_text('#!/bin/sh\n', encoding = 'utf-8')

    def upload_queue_job(self, model_name, runtime_name, file_list):
        self.uploads.append((model_name, runtime_name, file_list))

    def launch_queue_job(self, model_name, runtime_name, restart, batch):
        self.launches.append((model_name, runtime_name, restart, batch))


def _model(cluster, stagingpath):
    return SimpleNamespace(
        cluster = cluster,
        settings = SimpleNamespace(
            stagingpath = str(stagingpath),
            io_gather = 1,
            waitonlock = 0,
        ),
        private = SimpleNamespace(solution = '', runtimename = ''),
        miscellaneous = SimpleNamespace(name = 'testmodel'),
        verbose = SimpleNamespace(solution = False),
        qmu = SimpleNamespace(isdakota = False),
        toolkits = _FakeToolkits(),
        debug = SimpleNamespace(valgrind = False, gprof = False),
        transient = SimpleNamespace(isoceancoupling = False),
    )


def _fake_marshall(md, filename = None):
    Path(filename).write_bytes(b'bin')


def test_issmsettings_defaults_stagingpath_to_issm_execution(tmp_path, monkeypatch):
    monkeypatch.setattr(config, 'get_issm_dir', lambda: str(tmp_path))

    settings = issmsettings()

    assert settings.stagingpath == str(tmp_path / 'execution')
    assert 'stagingpath' in repr(settings)


def test_prepare_staging_directory_uses_executionpath_for_local_run(tmp_path, monkeypatch):
    executionpath = tmp_path / 'execution'
    stagingpath = tmp_path / 'staging'
    executionpath.mkdir()
    stagingpath.mkdir()
    monkeypatch.setattr(execute.tools.config, 'get_hostname', lambda: 'local')
    md = _model(_FakeCluster('LOCAL', executionpath), stagingpath)
    md.private.runtimename = 'run-1'

    run_directory = execute._prepare_staging_directory(md)

    assert run_directory == str(executionpath / 'run-1')
    assert Path(run_directory).is_dir()


def test_prepare_staging_directory_uses_and_cleans_remote_stagingpath(tmp_path, monkeypatch):
    executionpath = tmp_path / 'remote-execution'
    stagingpath = tmp_path / 'staging'
    stagingpath.mkdir()
    stale_directory = stagingpath / 'run-1'
    stale_directory.mkdir()
    (stale_directory / 'stale.bin').write_bytes(b'stale')
    monkeypatch.setattr(execute.tools.config, 'get_hostname', lambda: 'local')
    md = _model(_FakeCluster('remote', executionpath), stagingpath)
    md.private.runtimename = 'run-1'

    run_directory = execute._prepare_staging_directory(md)

    assert run_directory == str(stale_directory)
    assert list(stale_directory.iterdir()) == []


def test_prepare_staging_directory_rejects_missing_path(tmp_path, monkeypatch):
    monkeypatch.setattr(execute.tools.config, 'get_hostname', lambda: 'local')
    md = _model(_FakeCluster('remote', tmp_path / 'remote'), tmp_path / 'missing')
    md.private.runtimename = 'run-1'

    with pytest.raises(RuntimeError, match = 'Could not find staging directory'):
        execute._prepare_staging_directory(md)


def test_create_input_archive_contains_only_input_basenames(tmp_path):
    run_directory = tmp_path / 'run-1'
    run_directory.mkdir()
    inputs = [run_directory / 'model.bin', run_directory / 'model.queue']
    for path in inputs:
        path.write_text(path.name, encoding = 'utf-8')

    archive_path = _create_input_archive('run-1', [str(path) for path in inputs])

    assert archive_path == str(tmp_path / 'run-1.tar.gz')
    with tarfile.open(archive_path, 'r:gz') as archive:
        assert archive.getnames() == ['model.bin', 'model.queue']


def test_marshall_accepts_explicit_destination(tmp_path):
    md = SimpleNamespace(
        miscellaneous = SimpleNamespace(name = 'model'),
        verbose = SimpleNamespace(solution = False),
        model_class_names = lambda: [],
    )
    destination = tmp_path / 'run-1' / 'model.bin'
    destination.parent.mkdir()

    execute.marshall(md, destination)

    assert destination.is_file()
    assert destination.stat().st_size > 0


def test_solve_stages_local_inputs_without_uploading(tmp_path, monkeypatch):
    executionpath = tmp_path / 'execution'
    stagingpath = tmp_path / 'staging'
    executionpath.mkdir()
    stagingpath.mkdir()
    cluster = _FakeCluster('local', executionpath)
    md = _model(cluster, stagingpath)
    monkeypatch.setattr(execute.tools.config, 'get_hostname', lambda: 'local')
    monkeypatch.setattr(execute.tools.config, 'is_pc', lambda: False)
    monkeypatch.setattr(execute, 'marshall', _fake_marshall)

    result = execute.solve(
        md,
        'stressbalance',
        batch = True,
        check_consistency = False,
        runtime_name = False,
    )

    run_directory = executionpath / 'testmodel'
    assert result is md
    assert (run_directory / 'testmodel.bin').is_file()
    assert (run_directory / 'testmodel.toolkits').is_file()
    assert (run_directory / 'testmodel.queue').is_file()
    assert cluster.queue_options['file_prefix'] == str(run_directory / 'testmodel')
    assert cluster.uploads == []
    assert cluster.launches == [('testmodel', 'testmodel', None, True)]


def test_solve_stages_remote_inputs_and_uploads_full_paths(tmp_path, monkeypatch):
    stagingpath = tmp_path / 'staging'
    stagingpath.mkdir()
    cluster = _FakeCluster('remote', '/remote/execution')
    md = _model(cluster, stagingpath)
    monkeypatch.setattr(execute.tools.config, 'get_hostname', lambda: 'local')
    monkeypatch.setattr(execute.tools.config, 'is_pc', lambda: False)
    monkeypatch.setattr(execute, 'marshall', _fake_marshall)

    execute.solve(
        md,
        'stressbalance',
        batch = True,
        check_consistency = False,
        runtime_name = False,
    )

    run_directory = stagingpath / 'testmodel'
    assert len(cluster.uploads) == 1
    _, _, file_list = cluster.uploads[0]
    assert file_list == [
        str(run_directory / 'testmodel.bin'),
        str(run_directory / 'testmodel.toolkits'),
        str(run_directory / 'testmodel.queue'),
    ]


def test_restart_name_is_not_overwritten_by_runtime_name_option(tmp_path, monkeypatch):
    stagingpath = tmp_path / 'missing-staging'
    cluster = _FakeCluster('remote', '/remote/execution')
    md = _model(cluster, stagingpath)
    monkeypatch.setattr(execute.tools.config, 'get_hostname', lambda: 'local')

    def fail_marshall(*args, **kwargs):
        raise AssertionError('restart should not rewrite staged inputs')

    monkeypatch.setattr(execute, 'marshall', fail_marshall)

    execute.solve(
        md,
        'stressbalance',
        batch = True,
        check_consistency = False,
        restart = 'existing-run',
        runtime_name = True,
    )

    assert md.private.runtimename == 'existing-run'
    assert cluster.queue_options is None
    assert cluster.uploads == []
    assert cluster.launches == [('testmodel', 'existing-run', 'existing-run', True)]
