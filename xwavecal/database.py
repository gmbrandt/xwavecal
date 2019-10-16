import sqlite3
import os
from datetime import datetime

import logging
logger = logging.getLogger(__name__)


def add_data_to_db(db_path, db_info):
    db_exists = os.path.exists(db_path)
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    if not db_exists:
        logger.info('No existing database found at {0}. Creating a new one at {1}.'.format(db_path, db_path))
        # Create table
        c.execute('''CREATE TABLE caldata (id integer primary key autoincrement, 
                     type text, observation_date text,
                     date_created text, instrument text, instrument2 text, site_name text, 
                     fiber0 integer, fiber1 integer, fiber2 integer, 
                     is_bad integer, filepath text)''')
    match = query_db_for_match(c, db_info)
    if match is None:
        entries = tuple(db_info[key] for key in ['type', 'observation_date', 'date_created', 'instrument',
                                                 'instrument2',
                                                 'site_name', 'fiber0', 'fiber1', 'fiber2', 'is_bad',
                                                 'filepath'])
        c.execute('INSERT INTO caldata (type, observation_date, date_created, instrument, instrument2, site_name,'
                      ' fiber0, fiber1, fiber2, is_bad, filepath) '
                      'VALUES (?,?,?,?,?,?,?,?,?,?,?)', entries)
    else:
        c.execute('UPDATE caldata SET date_created=?, is_bad=?, filepath=? WHERE id=?',
                  (db_info['date_created'], db_info['is_bad'], db_info['filepath'], match))
    conn.commit()
    conn.close()


def format_db_info(data, time_format='%Y-%m-%dT%H:%M:%S.%f'):
    db_info = {'type': data.get_header_val('type'),
               'observation_date': datetime.strptime(data.get_header_val('observation_date'), time_format),
               'date_created': datetime.now(),
               'instrument': data.get_header_val('instrument'),
               'instrument2': data.get_header_val('instrument2'),
               'site_name': data.get_header_val('site_name'),
               'fiber0': int(data.fiber0_lit),
               'fiber1': int(data.fiber1_lit),
               'fiber2': int(data.fiber2_lit),
               'is_bad': 0, 'filepath': data.filepath}
    return db_info


def query_db_for_nearest(db_path, data, type_to_find, time_format='%Y-%m-%dT%H:%M:%S.%f', date=None):
    if not os.path.exists(db_path):
        logger.error('No existing database found at {0}. Cannot fetch calibration.'.format(db_path, db_path))
        return None
    db_info = format_db_info(data, time_format)
    if date is None:
        date = db_info['observation_date']
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    search_criteria = (type_to_find, db_info['instrument'], db_info['instrument2'], db_info['site_name'],
                       db_info['fiber0'], db_info['fiber1'], db_info['fiber2'],
                       date)
    c.execute('SELECT filepath FROM caldata WHERE is_bad=0 AND type=? AND instrument=? AND instrument2=?'
              'AND site_name=? AND fiber0=? AND fiber1=? AND fiber2=?'
              'ORDER BY ABS( strftime("%s", observation_date) - strftime("%s", ?) ) ASC',
               search_criteria)
    results = c.fetchall()
    conn.close()
    if len(results) == 0:
        return None
    else:
        return results[0][0]


def query_db_for_match(cursor, db_info):
    search_criteria = tuple(db_info[key] for key in ['type', 'observation_date',
                                                     'instrument', 'instrument2', 'site_name',
                                                     'fiber0', 'fiber1', 'fiber2'])
    cursor.execute('SELECT id FROM caldata WHERE type=? AND observation_date=?'
                   'AND instrument=? AND instrument2=? AND site_name=? AND fiber0=? AND fiber1=? AND fiber2=?',
                    search_criteria)
    results = cursor.fetchall()
    if len(results) == 0:
        return None
    else:
        return results[0][0]
