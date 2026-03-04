from climateset.download.utils import isolated_esgpull_context


def test_isolated_esgpull_context(tmp_path):
    with isolated_esgpull_context(tmp_path) as esg:
        assert esg is not None
        # the path should be tmp_path / .esgpull_jobs / <uuid>
        esg_path = esg.path
        assert esg_path.parent.name == ".esgpull_jobs"
        assert esg_path.parent.parent == tmp_path

        # It should exist during the context
        assert esg_path.exists()

    # After the context, it should be deleted
    assert not esg_path.exists()
