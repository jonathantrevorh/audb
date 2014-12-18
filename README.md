audb
====
## TODO
- iteratively build the index tree

## Running
Not final packaging, but running `./sqbuild` compiles the extension and run test commands in sqlite3.

## Running
 1. Install requirements as per below
 2. Run `sqlite3` to get the sqlite3 shell
 3. Load the binary by running `.load ./audb`
 4. Run the test sql with `.read test.sql` (it will fail because we've removed the music files from the public distribution)
 5. Optionally load your own music to try it out for yourself
 6. Run `.quit` to exit sqlite3

## Requirements
`sqlite3` and its development files are required to run and compile against sqlite3.
Installing on Ubuntu:
```
sudo apt-get install sqlite3 sqlite3-dev
```

## Building from source
Compile the file `audb.cc`. We've supplied a bash script which depends on `g++`, and passes the correct arguments for you. Just run `./build` from bash.
