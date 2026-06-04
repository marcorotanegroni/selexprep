# API reference

`selexprep`'s primary interface is the [CLI](cli.md). For programmatic use,
the two schemas below are the stable, documented public surface in v0.1: the
`LibraryReport` produced by primer inference, and the run provenance manifest.

See the [LibraryReport reference](library-report.md) for a narrative,
field-by-field walkthrough of inference outputs.

## LibraryReport

::: selexprep.library.report.LibraryReport
    options:
      show_root_heading: true
      show_source: false
      members_order: source

## SelexprepManifestV1

::: selexprep.manifest.SelexprepManifestV1
    options:
      show_root_heading: true
      show_source: false
      members_order: source
