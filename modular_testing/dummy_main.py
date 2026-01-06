import importlib
import logging

logger = logging.getLogger(__name__)

# Dummy job for testing
class DummyJob:
    def __init__(self, job_id, site_id, project_path, start_step, end_step, status="pending"):
        self.job_id = job_id
        self.site_id = site_id
        self.project_path = project_path
        self.start_step = start_step
        self.end_step = end_step
        self.status = status

def get_step_function(step_number: int):
    """Dynamically import and return the run function for the given step."""
    module_name = f"dummy_steps.step_{step_number:02d}_initialize" if step_number == 1 else f"dummy_steps.step_{step_number:02d}_process" if step_number == 2 else f"dummy_steps.step_{step_number:02d}_finalize"
    try:
        module = importlib.import_module(module_name)
        return module.run
    except ImportError as e:
        logger.error(f"Could not load module {module_name}: {e}")
        return None

def main():
    job = DummyJob(
        job_id=1,
        site_id="TEST_SITE",
        project_path="C:/tmp/test_project",
        start_step=1,
        end_step=3
    )

    for step in range(job.start_step, job.end_step + 1):
        logger.info(f"Running step {step}")
        step_func = get_step_function(step)
        if step_func:
            step_func(job)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
