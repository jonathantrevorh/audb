audb
====
## TODO
- vptree impl - done
- emailing prof to verify that we use sqlite instead what was proposed - done
- investigation on virtual table integration with sqlite - done (it's possible and the path we've chosen)
- implementing virtual table to connect to vp tree impl - in progress

## Requirements
`sqlite3` and its development files are required to run and compile against sqlite3.
Installing on Ubuntu:
```
sudo apt-get install sqlite3 sqlite3-dev
```

## Running
Not final packaging, but running `./sqbuild` compiles the extension and run test commands in sqlite3.
