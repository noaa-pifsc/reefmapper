```mermaid
erDiagram
    SITE_VISIT ||--o{ SFM_METADATA : contains
    SFM_BATCHES ||--o{ SFM_PROCESSING_JOBS : contains
    SFM_METADATA ||--o{ SFM_PROCESSING_JOBS : processes
    SFM_PROCESSING_JOBS ||--o{ SFM_PROCESSING_LOGS : generates

    SITE_VISIT {
        number SITEVISITID PK
        varchar2 SITE
        varchar2 ISLANDCODE
        varchar2 REGIONCODE
        number LATITUDE_N
        number LONGITUDE_E
        timestamp created_at
        timestamp updated_at
    }

    SFM_METADATA {
        number SFMMETAID PK
        number SITEVISITID FK
        varchar2 FILEPATH
        varchar2 MARKER_0M_NUMBERS
        varchar2 MARKER_5M_NUMBERS
        varchar2 MARKER_10M_NUMBERS
        number MARKER_0M_DEPTH_FT
        number MARKER_5M_DEPTH_FT
        number MARKER_10M_DEPTH_FT
        varchar2 SURVEY_TYPE
        number TOTAL_IMAGES
        timestamp created_at
        timestamp updated_at
    }

    SFM_BATCHES {
        number batch_id PK
        varchar2 batch_name
        varchar2 description
        timestamp created_at
        timestamp updated_at
        varchar2 status
    }

    SFM_PROCESSING_JOBS {
        number job_id PK
        number batch_id FK
        number sfmmetaid FK
        varchar2 project_path
        number start_step
        number end_step
        number quality
        varchar2 survey_year
        number priority
        varchar2 status
        varchar2 error_message
        timestamp created_at
        timestamp updated_at
        timestamp started_at
        timestamp completed_at
    }

    SFM_PROCESSING_LOGS {
        number log_id PK
        number job_id FK
        number step_num
        varchar2 step_name
        varchar2 status
        varchar2 message
        timestamp timestamp
    }
```

# Database Structure

## Tables

### SITE_VISIT
- Stores information about reef site visits
- Includes location, codes, and mission info

### SFM_METADATA
- Contains detailed metadata for each site visit
- Marker numbers, depths, survey type, total images, etc.

### SFM_BATCHES
- Groups multiple processing jobs together
- Tracks overall batch status and metadata

### SFM_PROCESSING_JOBS
- Main processing job information
- Tracks individual site processing status and progress
- Links to both batches and SFM_METADATA

### SFM_PROCESSING_LOGS
- Detailed processing logs
- Step-by-step status updates and messages
- Links to processing jobs

## Views

### V_SFM_ACTIVE_JOBS
- Shows currently running jobs
- Includes site and batch information
- Calculates running time

### V_SFM_JOB_STATS
- Job-level statistics
- Error and warning counts
- Processing duration

### V_SFM_STEP_STATS
- Statistics by processing step
- Average completion times
- Error counts by step

## Setup Instructions

1. Install required packages:
```bash
pip install cx-Oracle python-dotenv
```

2. Create the `.env` file:
```env
DB_TYPE=oracle
ORACLE_USER=your_username
ORACLE_PASSWORD=your_password
ORACLE_DSN=your_connection_string
ORACLE_MIN_CONNECTIONS=2
ORACLE_MAX_CONNECTIONS=5
ORACLE_CONNECTION_INCREMENT=1
```

3. Run the DDL scripts in order:
```sql
@01_create_tables.sql
@02_create_views.sql
```

## Database Relationships

- Each batch can contain multiple processing jobs
- Each SFM_METADATA record can have multiple processing jobs (for different years/surveys)
- Each site visit can have multiple metadata records
- Each job generates multiple log entries
- Jobs are linked to both batches and SFM_METADATA

## Status Values

### Batch Status
- pending
- processing
- completed
- failed

### Job Status
- pending
- running
- completed
- failed

### Log Status
- started
- running
- completed
- failed
- warning
- error
- info

## Notes

- All tables include created_at/updated_at timestamps
- Jobs track priority for processing order
- Marker numbers and depths are stored in SFM_METADATA
- Step numbers range from -1 to 7 (-1 for system messages)
- Quality values must be between 0 and 1
