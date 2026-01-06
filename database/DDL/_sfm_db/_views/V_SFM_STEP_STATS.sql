--------------------------------------------------------
--  DDL for View V_SFM_STEP_STATS
--------------------------------------------------------

  CREATE OR REPLACE EDITIONABLE VIEW "V_SFM_STEP_STATS" ("STEP_NAME", "JOB_COUNT", "AVG_DURATION_HOURS", "ERROR_COUNT") AS SELECT
    l_start.STEP_NAME,
    COUNT(DISTINCT l_start.JOB_ID) AS JOB_COUNT,
    ROUND(AVG((CAST(l_end.TIMESTAMP AS DATE) - CAST(l_start.TIMESTAMP AS DATE)) * 24), 2) AS AVG_DURATION_HOURS,
    COUNT(CASE 
        WHEN l_end.STATUS = 'FAILED' OR l_end.MESSAGE LIKE '%error%' 
        THEN 1 END
    ) AS ERROR_COUNT
FROM
    SFM_PROCESSING_LOGS l_start
JOIN
    SFM_PROCESSING_LOGS l_end
    ON l_start.JOB_ID = l_end.JOB_ID
    AND l_start.STEP_NAME = l_end.STEP_NAME
    AND l_start.STATUS = 'STARTED'
    AND (l_end.STATUS IN ('SUCCESS', 'FAILED') OR l_end.MESSAGE = 'Completed')
    AND l_end.TIMESTAMP > l_start.TIMESTAMP
WHERE
    l_start.TIMESTAMP IS NOT NULL
    AND l_end.TIMESTAMP IS NOT NULL
GROUP BY
    l_start.STEP_NAME
