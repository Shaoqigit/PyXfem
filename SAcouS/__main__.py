import sys
import argparse
import code
import readline

from SAcouS.interface.sol_setup import PyAcoustiXSetuper


class PyAcoustiXInteractiveConsole(code.InteractiveConsole):
    def __init__(self, locals=None, filename="<console>", histfile="~/.pyacoustix_history"):
        readline.parse_and_bind('tab: complete')
        readline.parse_and_bind('"\e[A": history-search-backward')
        readline.parse_and_bind('"\e[B": history-search-forward')
        super(PyAcoustiXInteractiveConsole, self).__init__(locals, filename)
        

    def interact(self, banner: str | None = None, exitmsg: str | None = None) -> None:
        try:
            sys.ps1
        except AttributeError:
            sys.ps1 = ">>> "
        try:
            sys.ps2
        except AttributeError:
            sys.ps2 = "... "
        cprt = 'Type "help", "copyright" or "license" for more information.\nAuthor: Shaoqiwu@outlook.com'
        if banner is None:
            self.write("PyAcoustiX 0.9.9\nbase on python %s on %s\n%s\n(%s)\n" %
                       (sys.version, sys.platform, cprt,
                        self.__class__.__name__))
        elif banner:
            self.write("%s\n" % str(banner))
        more = 0
        while 1:
            try:
                if more:
                    prompt = sys.ps2
                else:
                    prompt = sys.ps1
                try:
                    line = self.raw_input(prompt)
                except EOFError:
                    self.write("\n")
                    break
                else:
                    more = self.push(line)
            except KeyboardInterrupt:
                self.write("\nKeyboardInterrupt\n")
                self.resetbuffer()
                more = 0
        if exitmsg is None:
            self.write('now exiting %s...\n' % self.__class__.__name__)
        elif exitmsg != '':
            self.write('%s\n' % exitmsg)

    def raw_input(self, prompt=""):
        line = readline.get_line_buffer()
        if line:
            return input(prompt)
        else:
            return super().raw_input(prompt)

    

def main():
    argv = sys.argv[1:]
    if len(argv) == 0:
        vars = globals().copy()
        vars.update(locals())
        shell = PyAcoustiXInteractiveConsole(vars)
        print(f'Welcome to PyAcoustiX console!')
        shell.interact()
    else:
        parser = argparse.ArgumentParser(
            description="Acoustic analysis with PyAcoustiX",
            formatter_class=argparse.RawTextHelpFormatter,
        )

        parser.add_argument("--input_file", "--i", type=str, help="input file to be read from")
        args = parser.parse_args(argv)

        file_path = args.input_file
        sol_setuper = PyAcoustiXSetuper()
        sol_setuper.parse_input(file_path)
        sol_setuper.welcome()
        read_info = sol_setuper.sol_info
        for key, value in read_info.items():
            print(f'{key}: {value}')
        sol_setuper.exit()