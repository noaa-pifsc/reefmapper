import oracledb
from typing import Optional, List, Union
from dataclasses import dataclass
import os
import sqlite3
from datetime import datetime
from dotenv import load_dotenv
from pathlib import Path

@dataclass
class ProcessingJob:
    id: int
    sfmmetaid: int
    project_path: str
    site_id: str
    start_step: int
    end_step: int
    quality: float
    survey_year: str
    status: str
    priority: Optional[int] = None
    error_message: Optional[str] = None
    marker_pair1: Optional[str] = None
    marker_pair2: Optional[str] = None
    marker_pair3: Optional[str] = None
    marker_pair4: Optional[str] = None


class MyDatabaseManager:
    def __init__(self):
        # Load environment variables
        env_path = Path(__file__).parent.parent / '.env'
        load_dotenv(env_path)
        
        # Determine database type from environment
        self.use_oracle = os.getenv('DB_TYPE', 'sqlite').lower() == 'oracle'
        
        if self.use_oracle:
            # initialize Oracle client in thick mode for db credentials
            try:
                oracledb.init_oracle_client(lib_dir=r"C:\oracle\instantclient_21_19")
            except oracledb.OracleClientError as e:
                print(f"Failed to initialize Oracle Client: {e}")
           
            # Oracle connection setup from .env
            self.oracle_pool = oracledb.create_pool(
                user=os.getenv('ORACLE_USER'),
                password=os.getenv('ORACLE_PASSWORD'),
                dsn=os.getenv('ORACLE_DSN'),  # hostname:port/servicename format
                min=int(os.getenv('ORACLE_MIN_CONNECTIONS', '1')),
                max=int(os.getenv('ORACLE_MAX_CONNECTIONS', '3')),
                increment=int(os.getenv('ORACLE_CONNECTION_INCREMENT', '1'))
            )

        else:
            # SQLite setup from .env
            self.db_path = os.getenv('SQLITE_DB_PATH', 'database/reefmapper.db')
            os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
            self.conn = sqlite3.connect(self.db_path)
            self._init_sqlite_db()

    def connect(self):
        """Get appropriate database connection based on configuration."""
        if self.use_oracle:
            return self.oracle_pool.acquire()
        return self.conn

    def _release_connection(self, conn):
        """Release connection back to pool if using Oracle."""
        if self.use_oracle and conn:
            self.oracle_pool.release(conn)

    def close(self, conn=None):
        """Close or release a database connection."""
        if self.use_oracle:
            if conn:
                self._release_connection(conn)
            # Optionally, close the pool entirely (careful if you’re reusing)
            # self.oracle_pool.close()
        else:
            if self.conn:
                self.conn.close()
    
    def get_pending_jobs(self) -> List[ProcessingJob]:
        """Get all pending jobs ordered by priority and creation time."""
        conn = self.connect()
        try:
            cur = conn.cursor()
            if self.use_oracle:
                query = """
                    SELECT 
                    pj.job_id as id,
                    pj.sfmmetaid,
                    pj.project_path,
                    pj.site_id,
                    pj.start_step,
                    pj.end_step,
                    pj.quality,
                    pj.survey_year,
                    pj.status,
                    pj.priority,
                    pj.error_message,
                    sm.MARKER_0M_NUMBERS as marker_pair1,
                    sm.MARKER_10M_NUMBERS as marker_pair2,
                    sm.MARKER_15M_NUMBERS as marker_pair3,
                    sm.MARKER_20M_NUMBERS as marker_pair4
                    FROM SFM_PROCESSING_JOBS pj
                    JOIN sfm_metadata sm ON sm.sfmmetaid = pj.sfmmetaid
                    WHERE pj.status in ('pending', 'failed','completed')
                    AND pj.priority IS NOT NULL
                    and pj.end_step <>7
                    FETCH FIRST 1 ROWS ONLY
                """
                cur.execute(query)
            else:
                query = """
                    SELECT 
                    pj.job_id as id,
                    pj.sfmmetaid,
                    pj.project_path,
                    pj.site_id,
                    pj.start_step,
                    pj.end_step,
                    pj.quality,
                    pj.survey_year,
                    pj.status,
                    pj.priority,
                    pj.error_message,
                    sm.MARKER_0M_NUMBERS as marker_pair1,
                    sm.MARKER_10M_NUMBERS as marker_pair2,
                    sm.MARKER_15M_NUMBERS as marker_pair3,
                    sm.MARKER_20M_NUMBERS as marker_pair4
                    FROM SFM_PROCESSING_JOBS pj
                    JOIN sfm_metadata sm ON sm.sfmmetaid = pj.sfmmetaid
                    WHERE pj.status in ('pending', 'failed','completed')
                    AND pj.priority IS NOT NULL
                    and pj.end_step <>7
                    ORDER BY pj.priority ASC
                    FETCH FIRST 1 ROWS ONLY
                """
                cur.execute(query)

            jobs = []
            for row in cur.fetchall():
                jobs.append(ProcessingJob(*row))
            return jobs
        finally:
            if self.use_oracle:
                self._release_connection(conn)

    def log_step(self, job_id: int, step_num: int, status: str, message: str):
        """Log processing step with transaction support."""
        conn = self.connect()
        try:
            cur = conn.cursor()
            timestamp = datetime.now()
            step_name = f"Step {step_num}"
            if self.use_oracle:
                query = """
                    INSERT INTO SFM_PROCESSING_LOGS
                    (job_id, step_num, step_name, status, message, timestamp)
                    VALUES (:1, :2, :3, :4, :5, :6)
                """
                cur.execute(query, (job_id, step_num, step_name, status, message, timestamp))
            else:
                query = """
                    INSERT INTO SFM_PROCESSING_LOGS
                    (job_id, step_num, step_name, status, message, timestamp)
                    VALUES (?, ?, ?, ?, ?, ?)
                """
                cur.execute(query, (job_id, step_num, step_name, status, message, timestamp))
            conn.commit()
        except Exception as e:
            conn.rollback()
            raise
        finally:
            if self.use_oracle:
                self._release_connection(conn)

    def update_job_status(self, job_id: int, status: str, error_message: Optional[str] = None):
        """Update job status with transaction support."""
        conn = self.connect()
        try:
            cur = conn.cursor()
            timestamp = datetime.now()
            if self.use_oracle:
                query = """
                    UPDATE SFM_PROCESSING_JOBS
                    SET status = :1,
                        error_message = :2,
                        updated_at = :3
                    WHERE job_id = :4
                """
                cur.execute(query, (status, error_message, timestamp, job_id))
            else:
                query = """
                    UPDATE SFM_PROCESSING_JOBS
                    SET status = ?,
                        error_message = ?,
                        updated_at = ?
                    WHERE job_id = ?
                """
                cur.execute(query, (status, error_message, timestamp, job_id))
            conn.commit()
        except Exception as e:
            conn.rollback()
            raise
        finally:
            if self.use_oracle:
                self._release_connection(conn)

    def set_start_step(self, job_id: int, value: int):
        conn = self.connect()
        try:
            cur = conn.cursor()
            if self.use_oracle:
                query = "UPDATE sfm_processing_jobs SET start_step = :value WHERE job_id = :job_id"
                cur.execute(query, {"value": value, "job_id": job_id})
            else:
                query = "UPDATE sfm_processing_jobs SET start_step = ? WHERE job_id = ?"
                cur.execute(query, (value, job_id))
            conn.commit()
        except Exception as e:
            conn.rollback()
            raise
        finally:
            if self.use_oracle:
                self._release_connection(conn)

    def set_end_step(self, job_id: int, value: int):
        conn = self.connect()
        try:
            cur = conn.cursor()
            if self.use_oracle:
                query = "UPDATE sfm_processing_jobs SET end_step = :value WHERE job_id = :job_id"
                cur.execute(query, {"value": value, "job_id": job_id})
            else:
                query = "UPDATE sfm_processing_jobs SET end_step = ? WHERE job_id = ?"
                cur.execute(query, (value, job_id))
            conn.commit()
        except Exception as e:
            conn.rollback()
            raise
        finally:
            if self.use_oracle:
                self._release_connection(conn)

    def update_quality (self, job_id: int, value: float):
        conn = self.connect()
        try:
            cur = conn.cursor()
            if self.use_oracle:
                query = "UPDATE SFM_PROCESSING_JOBS SET quality = :1 WHERE job_id = :2"
                cur.execute(query, (value, job_id))
            else:
                query = "UPDATE SFM_PROCESSING_JOBS SET quality = ? WHERE job_id = ?"
                cur.execute(query, (value, job_id))
            conn.commit()
        except Exception as e:
            conn.rollback()
            raise
        finally:
            if self.use_oracle:
                self._release_connection(conn)

    def __del__(self):
        """Cleanup connections on object destruction."""
        if hasattr(self, 'oracle_pool'):
            self.oracle_pool.close()
        if hasattr(self, 'conn'):
            self.conn.close()

    def get_step_status(self, job_id: int, step_num: int) -> str:
        conn = self.connect()
        try:
            cur = conn.cursor()
            if self.use_oracle:
                query = """
                    SELECT status FROM SFM_PROCESSING_LOGS 
                    WHERE job_id = :1 AND step_num = :2 
                    ORDER BY timestamp DESC FETCH FIRST 1 ROWS ONLY
                """
                cur.execute(query, (job_id, step_num))
            else:
                query = """
                    SELECT status FROM SFM_PROCESSING_LOGS 
                    WHERE job_id = ? AND step_num = ? 
                    ORDER BY timestamp DESC LIMIT 1
                """
                cur.execute(query, (job_id, step_num))
            row = cur.fetchone()
            return row[0] if row else "not_started"
        finally:
            if self.use_oracle:
                self._release_connection(conn)
