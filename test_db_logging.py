import sqlite3
import oracledb
import datetime
import os
import time
import random
from dotenv import load_dotenv
from pathlib import Path

# Load environment variables
env_path = Path(__file__).parent / '.env'
load_dotenv(env_path)

# Determine database type from environment
USE_ORACLE = os.getenv('DB_TYPE', 'sqlite').lower() == 'oracle'

# Database configurations
SQLITE_PATH = os.getenv('SQLITE_DB_PATH', 'process_log_test.db')

# Oracle connection pool
oracle_pool = None
if USE_ORACLE:
    try:
        # Initialize Oracle client in thick mode for enterprise connections
        # This is required for connections with network encryption/security
        oracledb.init_oracle_client()
        print("Oracle thick mode initialized successfully")
    except Exception as e:
        print(f"Oracle thick mode initialization: {e}")
        print("Continuing with thin mode (may have connection limitations)")
    
    # Connection string format: hostname:port/servicename or hostname/servicename (if port 1521)
    # Examples: 
    #   - "localhost:1521/XEPDB1"
    #   - "db-server.company.com/ORCL"
    #   - "192.168.1.100:1522/PROD"
    oracle_pool = oracledb.create_pool(
        user=os.getenv('ORACLE_USER'),
        password=os.getenv('ORACLE_PASSWORD'),
        dsn=os.getenv('ORACLE_DSN'),  # hostname:port/servicename format
        min=int(os.getenv('ORACLE_MIN_CONNECTIONS', '1')),
        max=int(os.getenv('ORACLE_MAX_CONNECTIONS', '3')),
        increment=int(os.getenv('ORACLE_CONNECTION_INCREMENT', '1'))
    )

def get_connection():
    """Get database connection based on configuration."""
    if USE_ORACLE:
        return oracle_pool.acquire()
    return sqlite3.connect(SQLITE_PATH)

def release_connection(conn):
    """Release connection if using Oracle."""
    if USE_ORACLE and conn:
        oracle_pool.release(conn)

def init_db():
    """Check if database tables exist (Oracle) or create them (SQLite)."""
    if USE_ORACLE:
        conn = get_connection()
        try:
            c = conn.cursor()
            # Check if SFM_PROCESSING_LOGS table exists in Oracle
            c.execute("""
                SELECT table_name 
                FROM user_tables 
                WHERE table_name = 'SFM_PROCESSING_LOGS'
            """)
            table_exists = c.fetchone()
            if table_exists:
                print("Oracle: SFM_PROCESSING_LOGS table exists")
                
                # Check column definitions to understand size limits
                c.execute("""
                    SELECT column_name, data_type, data_length, data_precision
                    FROM user_tab_columns 
                    WHERE table_name = 'SFM_PROCESSING_LOGS'
                    AND column_name IN ('STEP_NAME', 'STATUS', 'MESSAGE')
                    ORDER BY column_name
                """)
                columns = c.fetchall()
                print("Column definitions:")
                for col in columns:
                    print(f"  {col[0]}: {col[1]}({col[2]})")
            else:
                print("Oracle: SFM_PROCESSING_LOGS table does not exist - please run DDL scripts")
                
            # Check if SFM_PROCESSING_JOBS table exists
            c.execute("""
                SELECT table_name 
                FROM user_tables 
                WHERE table_name = 'SFM_PROCESSING_JOBS'
            """)
            jobs_table_exists = c.fetchone()
            if jobs_table_exists:
                print("Oracle: SFM_PROCESSING_JOBS table exists")
                
                # Check ERROR_MESSAGE column size
                c.execute("""
                    SELECT column_name, data_type, data_length
                    FROM user_tab_columns 
                    WHERE table_name = 'SFM_PROCESSING_JOBS'
                    AND column_name = 'ERROR_MESSAGE'
                """)
                error_col = c.fetchone()
                if error_col:
                    print(f"  ERROR_MESSAGE: {error_col[1]}({error_col[2]})")
            else:
                print("Oracle: SFM_PROCESSING_JOBS table does not exist - please run DDL scripts")
        finally:
            release_connection(conn)
    else:
        conn = get_connection()
        try:
            c = conn.cursor()
            # Create SQLite tables for testing
            c.execute('''
                CREATE TABLE IF NOT EXISTS SFM_PROCESSING_JOBS (
                    JOB_ID INTEGER PRIMARY KEY AUTOINCREMENT,
                    BATCH_ID INTEGER,
                    SFMMETAID INTEGER,
                    PROJECT_PATH TEXT,
                    START_STEP INTEGER DEFAULT 1,
                    END_STEP INTEGER DEFAULT 7,
                    QUALITY REAL DEFAULT 0.5,
                    SURVEY_YEAR TEXT,
                    PRIORITY INTEGER DEFAULT 1,
                    STATUS TEXT DEFAULT 'pending',
                    ERROR_MESSAGE TEXT,
                    CREATED_AT TEXT,
                    UPDATED_AT TEXT,
                    STARTED_AT TEXT,
                    COMPLETED_AT TEXT
                )
            ''')
            c.execute('''
                CREATE TABLE IF NOT EXISTS SFM_PROCESSING_LOGS (
                    LOG_ID INTEGER PRIMARY KEY AUTOINCREMENT,
                    JOB_ID INTEGER,
                    STEP_NUM INTEGER,
                    STEP_NAME TEXT,
                    STATUS TEXT,
                    MESSAGE TEXT,
                    TIMESTAMP TEXT
                )
            ''')
            conn.commit()
            print("SQLite: Test tables created")
        finally:
            conn.close()

def log_step(job_id, step_num, step_name, status, message):
    """Log a processing step to the database."""
    conn = get_connection()
    try:
        c = conn.cursor()
        timestamp = datetime.datetime.now()
        
        if USE_ORACLE:
            query = """
                INSERT INTO SFM_PROCESSING_LOGS 
                (JOB_ID, STEP_NUM, STEP_NAME, STATUS, MESSAGE, TIMESTAMP)
                VALUES (:1, :2, :3, :4, :5, :6)
            """
        else:
            query = """
                INSERT INTO SFM_PROCESSING_LOGS 
                (JOB_ID, STEP_NUM, STEP_NAME, STATUS, MESSAGE, TIMESTAMP)
                VALUES (?, ?, ?, ?, ?, ?)
            """
            
        c.execute(query, (job_id, step_num, step_name, status, message, timestamp))
        conn.commit()
        print(f"[{timestamp.strftime('%H:%M:%S')}] Job {job_id} - Step {step_num}: {step_name} [{status}] {message}")
    finally:
        if USE_ORACLE:
            release_connection(conn)
        else:
            conn.close()

def update_job_status(job_id, status, error_message=None):
    """Update job status in the database."""
    conn = get_connection()
    try:
        c = conn.cursor()
        timestamp = datetime.datetime.now()
        
        # Truncate error message if it's too long for the Oracle column
        if error_message and len(error_message) > 20:  # Assuming 20 char limit based on error
            error_message = error_message[:17] + "..."
        
        if USE_ORACLE:
            query = """
                UPDATE SFM_PROCESSING_JOBS
                SET STATUS = :1, ERROR_MESSAGE = :2, UPDATED_AT = :3
                WHERE JOB_ID = :4
            """
        else:
            query = """
                UPDATE SFM_PROCESSING_JOBS
                SET STATUS = ?, ERROR_MESSAGE = ?, UPDATED_AT = ?
                WHERE JOB_ID = ?
            """
            
        c.execute(query, (status, error_message, timestamp, job_id))
        conn.commit()
        print(f"[{timestamp.strftime('%H:%M:%S')}] Job {job_id} status updated to: {status}")
    except Exception as e:
        print(f"Error updating job status: {e}")
    finally:
        if USE_ORACLE:
            release_connection(conn)
        else:
            conn.close()

def create_test_job():
    """Create a test job for processing."""
    conn = get_connection()
    try:
        c = conn.cursor()
        timestamp = datetime.datetime.now()
        
        if USE_ORACLE:
            query = """
                INSERT INTO SFM_PROCESSING_JOBS 
                (BATCH_ID, SFMMETAID, PROJECT_PATH, START_STEP, END_STEP, QUALITY, SURVEY_YEAR, STATUS, CREATED_AT)
                VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9)
            """
        else:
            query = """
                INSERT INTO SFM_PROCESSING_JOBS 
                (BATCH_ID, SFMMETAID, PROJECT_PATH, START_STEP, END_STEP, QUALITY, SURVEY_YEAR, STATUS, CREATED_AT)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """
            
        c.execute(query, (1, 999, '/test/path/to/images', 1, 7, 0.7, '2025', 'pending', timestamp))
        conn.commit()
        
        # Get the job ID
        if USE_ORACLE:
            c.execute("SELECT sfm_job_seq.CURRVAL FROM DUAL")
        else:
            c.execute("SELECT last_insert_rowid()")
        job_id = c.fetchone()[0]
        
        print(f"Created test job with ID: {job_id}")
        return job_id
    finally:
        if USE_ORACLE:
            release_connection(conn)
        else:
            conn.close()

def simulate_metashape_processing(job_id):
    """Simulate the actual Metashape processing steps with realistic timing."""
    
    # Processing steps based on the actual script - shortened names for Oracle compatibility
    steps = [
        (1, "Initialize", "Create project files and setup"),
        (2, "Align Photos", "Photo matching and alignment"),
        (3, "Filter Sparse", "Filter reconstruction uncertainty and projection accuracy"),
        (4, "Scale", "Detect markers and add scale bars"),
        (5, "Error Reduction", "Filter reprojection errors"),
        (6, "Dense Cloud", "Generate dense point cloud"),
        (7, "Export", "Export final products")
    ]
    
    print(f"\n=== Starting processing for Job {job_id} ===")
    update_job_status(job_id, 'running')
    log_step(job_id, 0, "Job Started", "started", f"Beginning batch processing for job {job_id}")
    
    for step_num, step_name, description in steps:
        print(f"\n--- Step {step_num}: {step_name} ---")
        
        # Log step start
        log_step(job_id, step_num, step_name, "running", f"Starting {description}")
        
        # Simulate processing time (different for each step)
        processing_times = {
            1: (2, 5),    # Initialize: 2-5 seconds
            2: (10, 20),  # Add and align photos: 10-20 seconds
            3: (5, 10),   # Sparse filtering: 5-10 seconds
            4: (8, 15),   # Scaling: 8-15 seconds
            5: (3, 7),    # Error reduction: 3-7 seconds
            6: (15, 30),  # Dense cloud: 15-30 seconds
            7: (12, 25)   # Export: 12-25 seconds
        }
        
        min_time, max_time = processing_times.get(step_num, (3, 8))
        processing_time = random.uniform(min_time, max_time)
        
        # Simulate sub-steps with progress updates
        sub_steps = max(2, int(processing_time / 3))
        for i in range(sub_steps):
            time.sleep(processing_time / sub_steps)
            progress = int((i + 1) / sub_steps * 100)
            if i < sub_steps - 1:  # Don't log 100% progress here
                log_step(job_id, step_num, step_name, "info", f"{description} - {progress}% complete")
        
        # Simulate occasional warnings or issues
        if random.random() < 0.3:  # 30% chance of warning
            warnings = [
                "Some images below quality threshold disabled",
                "Minor alignment issues detected but within acceptable range",
                "Marker detection required manual review",
                "Processing time longer than expected"
            ]
            warning_msg = random.choice(warnings)
            log_step(job_id, step_num, step_name, "warning", warning_msg)
        
        # Log step completion
        log_step(job_id, step_num, step_name, "completed", f"{step_name} finished successfully")
        print(f"✓ Step {step_num} completed in {processing_time:.1f} seconds")
    
    # Job completion
    update_job_status(job_id, 'completed')
    log_step(job_id, -1, "Job Completed", "completed", f"All processing steps completed successfully for job {job_id}")
    print(f"\n=== Job {job_id} completed successfully! ===")

def main():
    """Main test function."""
    print("=== SfM Database Logging Test ===")
    print(f"Database Type: {'Oracle' if USE_ORACLE else 'SQLite'}")
    
    if USE_ORACLE:
        print(f"Oracle DSN: {os.getenv('ORACLE_DSN')}")
        print(f"Oracle User: {os.getenv('ORACLE_USER')}")
        print("Note: Oracle thick mode is required for enterprise connections with encryption")
    
    init_db()
    
    if not USE_ORACLE:
        # Create a test job for SQLite
        job_id = create_test_job()
    else:
        # For Oracle, assume job already exists or create manually
        job_id = 1  # Use existing job ID or create one manually
        print(f"Using job ID: {job_id} (ensure this exists in Oracle)")
    
    # Simulate the processing
    try:
        simulate_metashape_processing(job_id)
    except Exception as e:
        print(f"Error during processing: {e}")
        update_job_status(job_id, 'failed', str(e))
        log_step(job_id, -1, "Job Failed", "error", f"Processing failed: {str(e)}")

if __name__ == "__main__":
    try:
        main()
    finally:
        # Clean up Oracle connection pool
        if oracle_pool:
            oracle_pool.close()
            print("\nOracle connection pool closed")
