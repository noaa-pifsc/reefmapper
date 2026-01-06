import cx_Oracle
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
    start_step: int
    end_step: int
    quality: float
    survey_year: str
    status: str
    error_message: Optional[str] = None

class DatabaseManager:
    def __init__(self):
        # Load environment variables
        env_path = Path(__file__).parent.parent / '.env'
        load_dotenv(env_path)
        
        # Determine database type from environment
        self.use_oracle = os.getenv('DB_TYPE', 'sqlite').lower() == 'oracle'
        
        if self.use_oracle:
            # Oracle connection setup from .env
            self.oracle_pool = cx_Oracle.SessionPool(
                user=os.getenv('ORACLE_USER'),
                password=os.getenv('ORACLE_PASSWORD'),
                dsn=os.getenv('ORACLE_DSN'),
                min=int(os.getenv('ORACLE_MIN_CONNECTIONS', '2')),
                max=int(os.getenv('ORACLE_MAX_CONNECTIONS', '5')),
                increment=int(os.getenv('ORACLE_CONNECTION_INCREMENT', '1')),
                encoding="UTF-8"
            )
        else:
            # SQLite setup from .env
            self.db_path = os.getenv('SQLITE_DB_PATH', 'database/reefmapper.db')
            os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
            self.conn = sqlite3.connect(self.db_path)
            self._init_sqlite_db()

    def _get_connection(self):
        """Get appropriate database connection based on configuration."""
        if self.use_oracle:
            return self.oracle_pool.acquire()
        return self.conn

    def _release_connection(self, conn):
        """Release connection back to pool if using Oracle."""
        if self.use_oracle and conn:
            self.oracle_pool.release(conn)

    def get_pending_jobs(self, batch_no: int) -> List[ProcessingJob]:
        """Get pending jobs for processing."""
        conn = self._get_connection()
        try:
            cur = conn.cursor()
            if self.use_oracle:
                query = """
                    SELECT job_id, sfmmetaid, project_path, start_step, end_step,
                           quality, survey_year, status, error_message
                    FROM SFM_PROCESSING_JOBS
                    WHERE batch_id = :1 AND status = 'pending'
                    ORDER BY priority DESC, created_at ASC
                """
                cur.execute(query, (batch_no,))
            else:
                query = """
                    SELECT job_id, sfmmetaid, project_path, start_step, end_step,
                           quality, survey_year, status, error_message
                    FROM SFM_PROCESSING_JOBS
                    WHERE batch_id = ? AND status = 'pending'
                    ORDER BY priority DESC, created_at ASC
                """
                cur.execute(query, (batch_no,))

            jobs = []
            for row in cur.fetchall():
                jobs.append(ProcessingJob(*row))
            return jobs
        finally:
            if self.use_oracle:
                self._release_connection(conn)

    def log_step(self, job_id: int, step_num: int, status: str, message: str):
        """Log processing step with transaction support."""
        conn = self._get_connection()
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
        conn = self._get_connection()
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

    def sync_to_oracle(self):
        """Sync local SQLite data to Oracle database."""
        if not self.use_oracle:
            return
        
        sqlite_conn = sqlite3.connect(self.db_path)
        oracle_conn = self._get_connection()
        
        try:
            # Sync jobs
            sqlite_cur = sqlite_conn.cursor()
            oracle_cur = oracle_conn.cursor()
            
            sqlite_cur.execute("SELECT * FROM processing_jobs WHERE synced = 0")
            for row in sqlite_cur.fetchall():
                # Insert into Oracle
                oracle_cur.execute("""
                    INSERT INTO reef_processing_jobs
                    VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11, :12, :13)
                """, row)
                
                # Mark as synced in SQLite
                sqlite_cur.execute("UPDATE processing_jobs SET synced = 1 WHERE id = ?", (row[0],))
            
            # Similar sync for logs table
            
            oracle_conn.commit()
            sqlite_conn.commit()
            
        except Exception as e:
            oracle_conn.rollback()
            sqlite_conn.rollback()
            raise
        finally:
            self._release_connection(oracle_conn)
            sqlite_conn.close()

    def __del__(self):
        """Cleanup connections on object destruction."""
        if hasattr(self, 'oracle_pool'):
            self.oracle_pool.close()
        if hasattr(self, 'conn'):
            self.conn.close()
