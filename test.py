import sqlite3
import time


def get_connection():
    conn = sqlite3.connect(":memory:")
    conn.enable_load_extension(True)
    conn.load_extension("./poc.so")
    conn.enable_load_extension(False)
    return conn

def get_cursor(conn):
    return conn.cursor()

def create_table(cursor):
    cursor.execute("create virtual table songs using audb_tree")

def insert_song(cursor, sid, path, song, artist):
    path = "\"" + path + "\""
    cursor.execute("insert into songs values(?, ?, ?, ?)", (sid, path, song, artist,))

def select_where(cursor, path):
    path = "\"" + path + "\""
    rows = cursor.execute("select * from songs where path = ?", (path,))
    return rows

def print_data(data):
    for row in data:
        print row

def close_connection(conn):
    conn.close()

def build_and_test(music, queries):
    # create connection, get cursor, set up table
    conn = get_connection()
    cursor = get_cursor(conn)
    create_table(cursor)

    # insert songs
    sid = 1
    for artist, songs in music.iteritems():
        for song in songs:
            insert_song(cursor, sid, song['path'], song['name'], artist)
            sid += 1


    timings = {}
    # search
    for query in queries:
        start = time.clock()
        data = select_where(cursor, query)
        end = time.clock()
        timings[query] = end - start

    print "average query time: {}".format(average(timings.values()))
    print "individual queries:"
    for query, timing in timings.iteritems():
        print "{} {}".format(timing, query)

    # clean up
    close_connection(conn)

def average(items):
    return sum(items) / len(items)

build_and_test({
    'Steve Aoki': [{
        'path': './Music/Aoki - Boneless.mp3',
        'name': 'Boneless'
    }],
    'AWOLNATION': [{
        'path': './Music/AWOLNATION - Sail.mp3',
        'name': 'Sail'
    }],
    'SKRILLEX': [{
        'path': './Music/Bangarang.mp3',
        'name': 'Bangarang'
    }],
    'Red Hot Chili Peppers': [{
        'path': './Music/RHCP - Look Around.mp3',
        'name': 'Look Around'
    }]
}, ['./Music/Aoki2.wav', './Music/RHCP3.wav'])
