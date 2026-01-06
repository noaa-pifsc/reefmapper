-- Create views for easier reporting and analysis

CREATE OR REPLACE VIEW V_SFM_ACTIVE_JOBS AS
SELECT 
    j.job_id,
    b.batch_name,
    m.FILEPATH,
    v.SITE,
    v.LATITUDE_N,
    v.LONGITUDE_E,
    j.start_step,
    j.end_step,
    j.status,
    j.created_at,
    j.started_at,
    ROUND((SYSDATE - j.started_at) * 24, 2) as hours_running
FROM 
    SFM_PROCESSING_JOBS j
    JOIN SFM_BATCHES b ON j.batch_id = b.batch_id
    JOIN SFM_METADATA m ON j.sfmmetaid = m.sfmmetaid
    JOIN SITE_VISIT v ON m.sitevisitid = v.sitevisitid
WHERE 
    j.status = 'running';

CREATE OR REPLACE VIEW V_SFM_JOB_STATS AS
SELECT 
    j.job_id,
    v.SITE,
    j.survey_year,
    j.status,
    j.started_at,
    j.completed_at,
    COUNT(l.log_id) as total_log_entries,
    SUM(CASE WHEN l.status = 'error' THEN 1 ELSE 0 END) as error_count,
    SUM(CASE WHEN l.status = 'warning' THEN 1 ELSE 0 END) as warning_count,
    ROUND((j.completed_at - j.started_at) * 24, 2) as total_hours
FROM 
    SFM_PROCESSING_JOBS j
    JOIN SFM_METADATA m ON j.sfmmetaid = m.sfmmetaid
    JOIN SITE_VISIT v ON m.sitevisitid = v.sitevisitid
    LEFT JOIN SFM_PROCESSING_LOGS l ON j.job_id = l.job_id
GROUP BY 
    j.job_id,
    v.SITE,
    j.survey_year,
    j.status,
    j.started_at,
    j.completed_at;

-- View for step completion statistics
CREATE OR REPLACE VIEW V_SFM_STEP_STATS AS
SELECT 
    l.step_num,
    l.step_name,
    COUNT(*) as total_executions,
    ROUND(AVG(CASE 
        WHEN l.status = 'completed' 
        THEN (LEAD(l.timestamp) OVER (PARTITION BY l.job_id, l.step_num ORDER BY l.timestamp) - l.timestamp) * 24 * 60
    END), 2) as avg_minutes_to_complete,
    SUM(CASE WHEN l.status = 'error' THEN 1 ELSE 0 END) as error_count
FROM 
    SFM_PROCESSING_LOGS l
GROUP BY 
    l.step_num,
    l.step_name
ORDER BY 
    l.step_num;
