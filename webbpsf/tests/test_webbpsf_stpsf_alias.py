def test_import_alias():
    import webbpsf
    import stpsf
    assert webbpsf.NIRCam == stpsf.NIRCam  # and repeat for the others too