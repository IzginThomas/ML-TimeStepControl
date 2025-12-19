# How to use

To run numerical experiments in parallel on a server to gather data,
you can execute something like

```bash
matlab -nodisplay -r "write_DATAMPRK22; exit"
```

and similar for the other `write_DATA*` files. Please check the number
of processes set in `parpool` in these files and adapt them to your system.
