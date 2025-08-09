--------------------------------------------------------
--  Constraints for Table SFM_PROCESSING_STEPS
--------------------------------------------------------

  ALTER TABLE "SFM_PROCESSING_STEPS" MODIFY ("DESCRIPTION" NOT NULL ENABLE)
  ALTER TABLE "SFM_PROCESSING_STEPS" ADD PRIMARY KEY ("STEP_NUM") USING INDEX  ENABLE
