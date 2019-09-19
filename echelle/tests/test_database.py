import tempfile
import os
from datetime import datetime, timedelta
import sqlite3

import echelle.database as db
from echelle.tests.utils import FakeImage


def test_format_db_info():
    data = FakeImage()
    fmt = '%Y-%m-%dT%H:%M:%S.%f'
    db_info = db.format_db_info(data, fmt)
    obstype, dateobs, datenow, inst, site, fib0, fib1, fib2, is_bad, path = db_info[0]
    assert all([obstype == data.get_header_val('type'),
                dateobs == datetime.strptime(data.get_header_val('observation_date'), fmt),
                datenow - datetime.now() <= timedelta(seconds=20),
                inst == data.get_header_val('instrument'),
                site == data.get_header_val('site_name'),
                fib0 == data.fiber0_lit, fib1 == data.fiber1_lit, fib2 == data.fiber2_lit,
                is_bad == 0, path == data.filepath])


def test_add_caldata_to_db():
    data = FakeImage()
    fmt = '%Y-%m-%dT%H:%M:%S.%f'
    db_info = db.format_db_info([data, data], fmt)
    with tempfile.TemporaryDirectory() as temp_directory:
        db_path = os.path.join(temp_directory, 'test.db')
        db.add_data_to_db(db_path, db_info)
        db.add_data_to_db(db_path, db_info[0])
        # expect 3 identical entries in the database.
        conn = sqlite3.connect(db_path)
        c = conn.cursor()
        c.execute('SELECT * FROM caldata')
        assert len(c.fetchall()) == 3
        conn.close()


def test_query_database():
    data1, data2, data3 = FakeImage(), FakeImage(), FakeImage()
    reference_data = FakeImage()
    data3.set_header_val('type', 'wrong')
    data3.set_header_val('observation_date', '2019-04-15T12:56:44.466')
    data2.set_header_val('observation_date', '2019-04-13T12:56:44.466')
    data2.filepath = 'some/path'
    data1.set_header_val('observation_date', '2019-04-11T12:56:44.466')
    reference_data.set_header_val('observation_date', '2019-04-16T12:56:44.466')
    fmt = '%Y-%m-%dT%H:%M:%S.%f'
    db_info = db.format_db_info([data1, data2, data3], fmt)
    with tempfile.TemporaryDirectory() as temp_directory:
        db_path = os.path.join(temp_directory, 'test.db')
        db.add_data_to_db(db_path, db_info)
        filepath = db.query_db_for_nearest(db_path, reference_data, 'lampflat', '%Y-%m-%dT%H:%M:%S.%f')
        assert filepath == data2.filepath
        assert db.query_db_for_nearest(db_path, reference_data, 'no_type', '%Y-%m-%dT%H:%M:%S.%f') is None


def test_query_database_returns_none():
    assert db.query_db_for_nearest('non/existent/path', None, 'lampflat', '%Y-%m-%dT%H:%M:%S.%f') is None
    assert db.query_db_for_nearest('non/existent/path', None, 'lampflat', '%Y-%m-%dT%H:%M:%S.%f') is None
