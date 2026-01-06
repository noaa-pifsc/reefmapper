--------------------------------------------------------
--  Constraints for Table SFM_PROCESSING_LOGS
--------------------------------------------------------

  ALTER TABLE "SFM_PROCESSING_LOGS" ADD CONSTRAINT "SFM_PROCESSING_LOGS_PK" PRIMARY KEY ("LOG_ID") USING INDEX  ENABLE
  ALTER TABLE "SFM_PROCESSING_LOGS" MODIFY ("LOG_ID" NOT NULL ENABLE)
