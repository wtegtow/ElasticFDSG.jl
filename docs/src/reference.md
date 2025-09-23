## API Reference

```@contents
Pages = ["reference.md"]
```

```@index
Pages = ["reference.md"]
```

```
@autodocs
Modules = [ElasticFDSG]  
```
```
@autodocs
Modules = [ElasticFDSG.dim2, ElasticFDSG.dim3]
filter = m -> occursin("runsim", string(m)) || occursin("solve!", string(m))
```