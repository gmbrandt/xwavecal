import tempfile
import os
from datetime import datetime, timedelta
import sqlite3
import mock

import xwavecal.database as db
from xwavecal.tests.utils import FakeImage


def test_format_db_info():
    data = FakeImage()
    fmt = '%Y-%m-%dT%H:%M:%S.%f'
    db_info = db.format_db_info(data, fmt)
    fib0, fib1, fib2 = db_info['fiber0'], db_info['fiber1'], db_info['fiber2']
    assert all([db_info['type'] == data.get_header_val('type'),
                db_info['observation_date'] == datetime.strptime(data.get_header_val('observation_date'), fmt),
                db_info['date_created'] - datetime.now() <= timedelta(seconds=20),
                db_info['instrument'] == data.get_header_val('instrument'),
                db_info['instrument2'] == data.get_header_val('instrument2'),
                db_info['site_name'] == data.get_header_val('site_name'),
                fib0 == data.fiber0_lit, fib1 == data.fiber1_lit, fib2 == data.fiber2_lit,
                db_info['is_bad'] == 0, db_info['filepath'] == data.filepath])


@mock.patch('xwavecal.database.query_db_for_match', return_value=None)
def test_add_to_db(mock_match):
    data = FakeImage()
    fmt = '%Y-%m-%dT%H:%M:%S.%f'
    db_info = db.format_db_info(data, fmt)
    with tempfile.TemporaryDirectory() as temp_directory:
        db_path = os.path.join(temp_directory, 'test.db')
        db.add_data_to_db(db_path, db_info)
        db.add_data_to_db(db_path, db_info)
        db.add_data_to_db(db_path, db_info)
        # expect 3 identical entries in the database.
        conn = sqlite3.connect(db_path)
        c = conn.cursor()
        c.execute('SELECT * FROM caldata')
        assert len(c.fetchall()) == 3
        conn.close()


def test_update_db():
    data = FakeImage()
    fmt = '%Y-%m-%dT%H:%M:%S.%f'
    db_info = db.format_db_info(data, fmt)
    with tempfile.TemporaryDirectory() as temp_directory:
        db_path = os.path.join(temp_directory, 'test.db')
        db.add_data_to_db(db_path, db_info)
        db.add_data_to_db(db_path, db_info)
        db.add_data_to_db(db_path, db_info)
        # expect the 3 identical entries in the database to be culled into one.
        conn = sqlite3.connect(db_path)
        c = conn.cursor()
        c.execute('SELECT * FROM caldata')
        assert len(c.fetchall()) == 1
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
    with tempfile.TemporaryDirectory() as temp_directory:
        db_path = os.path.join(temp_directory, 'test.db')
        for data in [data1, data2, data3]:
            db_info = db.format_db_info(data, fmt)
            db.add_data_to_db(db_path, db_info)
        filepath = db.query_db_for_nearest(db_path, reference_data, 'lampflat', fmt)
        assert filepath == data2.filepath
        assert db.query_db_for_nearest(db_path, reference_data, 'no_type', fmt) is None


def test_query_database_returns_none():
    assert db.query_db_for_nearest('non/existent/path', None, 'lampflat', '%Y-%m-%dT%H:%M:%S.%f') is None
    assert db.query_db_for_nearest('non/existent/path', None, 'lampflat', '%Y-%m-%dT%H:%M:%S.%f') is None
