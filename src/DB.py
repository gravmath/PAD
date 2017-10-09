from   PAD_settings import *
import sqlite3

instrument_conn   = sqlite3.connect(PAD_db_dir+'instrument_data_db.sqlite')
instrument_cursor = instrument_conn.cursor()
