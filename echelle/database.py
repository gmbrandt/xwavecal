import sqlite3
import os

import logging as logger
from datetime import datetime


def add_data_to_db(db_path, db_info):
    # note if the TABLE query can be a variable, this can be easily extended to any table in the database.
    if not isinstance(db_info, list):
        db_info = [db_info]
    db_exists = os.path.exists(db_path)
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    if not db_exists:
        logger.info('No existing database found at {0}. Creating a new one at {1}.'.format(db_path, db_path))
        # Create table
        c.execute('''CREATE TABLE caldata (type text, dateobs text, datecreated text, instrument text,
         site text, fiber0 integer, fiber1 integer, fiber2 integer,
          is_bad integer, filepath text)''')

    c.executemany('INSERT INTO caldata VALUES (?,?,?,?,?,?,?,?,?,?)', db_info)
    conn.commit()
    conn.close()


def format_db_info(all_data, time_format='%Y-%m-%dT%H:%M:%S.%f'):
    if not isinstance(all_data, list):
        all_data = [all_data]
    db_info = []
    for data in all_data:
        db_info.append((data.get_header_val('type'), datetime.strptime(data.get_header_val('observation_date'), time_format),
                        datetime.now(), data.get_header_val('instrument'),
                        data.get_header_val('site_name'),
                        int(data.fiber0_lit), int(data.fiber1_lit), int(data.fiber2_lit),
                        0, data.filepath))
    return db_info


def query_db_for_nearest(db_path, data, type_to_find, time_format='%Y-%m-%dT%H:%M:%S.%f', date=None):
    if not os.path.exists(db_path):
        logger.error('No existing database found at {0}. Cannot fetch calibration.'.format(db_path, db_path))
        return None
    if date is None:
        date = datetime.strptime(data.get_header_val('observation_date'), time_format)
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    search_criteria = (0, type_to_find, data.get_header_val('instrument'),
                       data.get_header_val('site_name'),
                       int(data.fiber0_lit), int(data.fiber1_lit), int(data.fiber2_lit),
                       date)
    c.execute('SELECT filepath FROM caldata WHERE is_bad=? AND type=? AND instrument=?'
              'AND site=? AND fiber0=? AND fiber1=? AND fiber2=?'
              'ORDER BY ABS( strftime("%s", dateobs) - strftime("%s", ?) ) ASC',
               search_criteria)
    results = c.fetchall()
    conn.close()
    if len(results) == 0:
        return None
    else:
        return results[0][0]