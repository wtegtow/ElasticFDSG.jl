## API Reference

```@contents
Pages = ["reference.md"]
```

```@index
Pages = ["reference.md"]
```

```@autodocs
Modules = [ElasticFDSG]
filter = f -> occursin("load_results", string(f)) ||
             occursin("print_h5_tree", string(f)) ||
             occursin("config_template", string(f))
```

```@autodocs
Modules = [ElasticFDSG.dim2, ElasticFDSG.dim3]
filter = f -> occursin("runsim", string(f))
```