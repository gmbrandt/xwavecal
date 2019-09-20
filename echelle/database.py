import sqlite3
import os

import logging as logger
from datetime import datetime


def add_data_to_db(db_path, db_info):
    db_exists = os.path.exists(db_path)
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    if not db_exists:
        logger.info('No existing database found at {0}. Creating a new one at {1}.'.format(db_path, db_path))
        # Create table
        c.execute('''CREATE TABLE caldata (id integer primary key autoincrement, type text, dateobs text,
                     datecreated text, instrument text, site text, 
                     fiber0 integer, fiber1 integer, fiber2 integer, 
                     is_bad integer, filepath text)''')
    match = query_db_for_match(c, db_info)
    if match is None:
        c.execute('INSERT INTO caldata (type, dateobs, datecreated, instrument, site,'
                      ' fiber0, fiber1, fiber2, is_bad, filepath) '
                      'VALUES (?,?,?,?,?,?,?,?,?,?)', db_info)
    else:
        c.execute('UPDATE caldata SET datecreated=?, is_bad=?, filepath=? WHERE id=?',
                  (db_info[2], db_info[8], db_info[9], match))
    conn.commit()
    conn.close()


def format_db_info(data, time_format='%Y-%m-%dT%H:%M:%S.%f'):
    db_info = ((data.get_header_val('type'), datetime.strptime(data.get_header_val('observation_date'), time_format),
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


def query_db_for_match(cursor, db_info):
    search_criteria = tuple(db_info[i] for i in [0, 1, 3, 4, 5, 6, 7])
    cursor.execute('SELECT id FROM caldata WHERE type=? AND dateobs=?'
                   'AND instrument=? AND site=? AND fiber0=? AND fiber1=? AND fiber2=?',
                    search_criteria)
    results = cursor.fetchall()
    if len(results) == 0:
        return None
    else:
        return results[0][0]
