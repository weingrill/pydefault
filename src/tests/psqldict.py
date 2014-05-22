'''
Created on May 22, 2014

@author: jwe

code to test the dictionary record with psql
'''
import psycopg2
import psycopg2.extras

database = 'corot'
user = 'sro'
host = 'pina'
conn = psycopg2.connect(database=database, user=user, host=host)
dict_cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
dict_cur.execute("SELECT * FROM hydra")
rec = dict_cur.fetchone()
print rec['object'], rec['image'], rec['vhelio']

conn.close()