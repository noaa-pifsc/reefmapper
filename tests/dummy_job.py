class DummyJob:
    def __init__(self):
        self.job_id = 1
        self.site_id = "TESTSITE"
        self.project_path = "/fake/path/project.psx"
        self.start_step = 1
        self.end_step = 3
        self.status = "pending"

    def __repr__(self):
        return f"<DummyJob id={self.job_id} site={self.site_id}>"
