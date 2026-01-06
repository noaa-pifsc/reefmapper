--------------------------------------------------------
--  Constraints for Table SFM_PROCESSING_JOBS
--------------------------------------------------------

  ALTER TABLE "SFM_PROCESSING_JOBS" ADD CONSTRAINT "SFM_PROCESSING_JOBS_PK" PRIMARY KEY ("JOB_ID") USING INDEX  ENABLE
  ALTER TABLE "SFM_PROCESSING_JOBS" MODIFY ("JOB_ID" NOT NULL ENABLE)
