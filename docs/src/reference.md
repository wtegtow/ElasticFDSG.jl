## API Reference

```@contents
Pages = ["reference.md"]
```

```@index
Pages = ["reference.md"]
```

```@autodocs
Modules = [ElasticFDSG, ElasticFDSG.dim2, ElasticFDSG.dim3]
Filter = f -> (
    string(f) in (
        "ElasticFDSG.config_template",
        "ElasticFDSG.load_results",
        "ElasticFDSG.print_h5_tree",
        "ElasticFDSG.dim2.runsim",
        "ElasticFDSG.dim3.runsim"
    )
)
```
