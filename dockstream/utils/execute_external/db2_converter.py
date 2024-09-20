from dockstream.utils.execute_external.execute import ExecutorBase


class DB2_converter_Executor(ExecutorBase):
    """For the execution of "DB2_converter_Executor" binaries."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(prefix_execution=prefix_execution, binary_location=binary_location)

    def execute(self, command: str, arguments: list, check=True, location=None):
        return super().execute(command=command,
                               arguments=arguments,
                               check=check,
                               location=location)

    def is_available(self):
        try:
            # self.execute(command="echo 1", arguments=[])
            # import os
            # print(os.environ["db2_converter_SH"])
            # print("db2_converter_SH")
            return True
        except Exception as e:
            print(e)
            return False
