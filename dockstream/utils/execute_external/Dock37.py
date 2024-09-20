from dockstream.utils.execute_external.execute import ExecutorBase


class Dock37Executor(ExecutorBase):
    """For the execution of "Dock37Executor" binaries."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(prefix_execution=prefix_execution, binary_location=binary_location)

    def execute(self, command: str, arguments: list, check=True, location=None):
        return super().execute(command=command,
                               arguments=arguments,
                               check=check,
                               location=location)

    def is_available(self):
        try:
            return True
        except Exception as e:
            return False
