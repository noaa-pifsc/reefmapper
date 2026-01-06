import sqlite3
from dataclasses import dataclass
from typing import List, Optional
from datetime import datetime
import os

@dataclass
class ProcessingJob:
    project_path: str
    site_id: str
    batch_code: int
    start_step: int
    end_step: int
    marker_pair1: str
    marker_pair2: str
    marker_pair3: str
    marker_pair4: str
    quality: float
    survey_year: str
    status: str
    log_notes: str = ""
    id: Optional[int] = None
    last_updated: Optional[str] = None

class DatabaseManager:
    def __init__(self, db_path: str = "reefmapper.db"):
        """Initialize database connection and create tables if they don't exist."""
        self.db_path = db_path
        self._create_tables()

    def _create_tables(self):
        """Create the necessary database tables if they don't exist."""
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()

        # Create processing jobs table
        c.execute('''
            CREATE TABLE IF NOT EXISTS processing_jobs (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                project_path TEXT NOT NULL,
                site_id TEXT NOT NULL,
                batch_code INTEGER NOT NULL,
                start_step INTEGER NOT NULL,
                end_step INTEGER NOT NULL,
                marker_pair1 TEXT,
                marker_pair2 TEXT,
                marker_pair3 TEXT,
                marker_pair4 TEXT,
                quality REAL NOT NULL,
                survey_year TEXT NOT NULL,
                status TEXT NOT NULL,
                log_notes TEXT,
                last_updated TEXT
            )
        ''')

        # Create processing logs table for detailed step logging
        c.execute('''
            CREATE TABLE IF NOT EXISTS processing_logs (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                job_id INTEGER NOT NULL,
                step INTEGER NOT NULL,
                status TEXT NOT NULL,
                message TEXT,
                timestamp TEXT NOT NULL,
                FOREIGN KEY (job_id) REFERENCES processing_jobs (id)
            )
        ''')

        conn.commit()
        conn.close()

    def add_job(self, job: ProcessingJob) -> int:
        """Add a new processing job to the database."""
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        
        c.execute('''
            INSERT INTO processing_jobs (
                project_path, site_id, batch_code, start_step, end_step,
                marker_pair1, marker_pair2, marker_pair3, marker_pair4,
                quality, survey_year, status, log_notes, last_updated
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            job.project_path, job.site_id, job.batch_code, job.start_step,
            job.end_step, job.marker_pair1, job.marker_pair2, job.marker_pair3,
            job.marker_pair4, job.quality, job.survey_year, job.status,
            job.log_notes, datetime.now().isoformat()
        ))
        
        job_id = c.lastrowid
        conn.commit()
        conn.close()
        return job_id

    def get_pending_jobs(self, batch_code: int) -> List[ProcessingJob]:
        """Get all pending jobs for a specific batch code."""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        
        c.execute('''
            SELECT * FROM processing_jobs
            WHERE batch_code = ? AND status = 'pending'
            ORDER BY id ASC
        ''', (batch_code,))
        
        rows = c.fetchall()
        conn.close()
        
        return [ProcessingJob(**dict(row)) for row in rows]

    def update_job_status(self, job_id: int, status: str, notes: str = None):
        """Update the status of a processing job."""
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        
        update_values = {
            'status': status,
            'last_updated': datetime.now().isoformat()
        }
        if notes:
            update_values['log_notes'] = notes
        
        set_clause = ', '.join(f'{k} = ?' for k in update_values.keys())
        query = f'UPDATE processing_jobs SET {set_clause} WHERE id = ?'
        
        c.execute(query, (*update_values.values(), job_id))
        conn.commit()
        conn.close()

    def log_step(self, job_id: int, step: int, status: str, message: str = None):
        """Log a processing step for a job."""
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        
        c.execute('''
            INSERT INTO processing_logs (job_id, step, status, message, timestamp)
            VALUES (?, ?, ?, ?, ?)
        ''', (job_id, step, status, message, datetime.now().isoformat()))
        
        conn.commit()
        conn.close()

    def import_from_csv(self, csv_path: str):
        """Import processing jobs from a CSV file."""
        import csv
        
        with open(csv_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                job = ProcessingJob(
                    project_path=row['Project File Path'],
                    site_id=row['Site ID'],
                    batch_code=int(row['Code']),
                    start_step=int(row['Start']),
                    end_step=int(row['End']),
                    marker_pair1=row['Marker Pair 1'],
                    marker_pair2=row['Marker Pair 2'],
                    marker_pair3=row['Marker Pair 3'],
                    marker_pair4=row['Marker Pair 4'],
                    quality=float(row['Quality']),
                    survey_year=row['Year'],
                    status='pending',
                    log_notes=row.get('log notes', '')
                )
                self.add_job(job)
