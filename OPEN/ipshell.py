# -*- coding: utf-8 -*-
#
# Copyright (c) 2013 João Faria
# This file is part of OPEN which is licensed under the MIT license.
# You should have received a copy of the license along with OPEN. See LICENSE.
#

"""
Provides a custom embedded IPython shell with its own magics.
"""

# standard library imports
import sys
import warnings

# ipython imports
# We need to use nested to support python 2.6, once we move to >=2.7, we can
# use the with keyword's new builtin support for nested managers
try:
    from contextlib import nested
except:
    from IPython.utils.nested_context import nested
from IPython.config.loader import Config
from IPython.core import ultratb
from IPython.core.magic import Magics, magics_class, line_magic
try: 
    # if IPython > 1.0 this works
    from IPython.terminal.interactiveshell import TerminalInteractiveShell
    from IPython.terminal.ipapp import load_default_config
except ImportError:
    # for IPython 0.13 we had to do this instead
    from IPython.frontend.terminal.interactiveshell import TerminalInteractiveShell
    from IPython.frontend.terminal.ipapp import load_default_config
from IPython.utils.traitlets import Bool, CBool, Unicode
from IPython.utils.io import ask_yes_no

# intra-package imports
from .commandsOPEN import EmbeddedMagics


class EmbedShell(TerminalInteractiveShell):

    dummy_mode = Bool(False)
    exit_msg = Unicode('')
    embedded = CBool(True)
    embedded_active = CBool(True)
    # Like the base class display_banner is not configurable, but here it
    # is True by default.
    display_banner = CBool(True)

    def __init__(self, config=None, ipython_dir=None, user_ns=None,
                 user_module=None, custom_exceptions=((),None),
                 usage=None, banner1=None, banner2=None,
                 display_banner=True, exit_msg=u'', user_global_ns=None):
    
        if user_global_ns is not None:
            warnings.warn("user_global_ns has been replaced by user_module. The\
                           parameter will be ignored.", DeprecationWarning)

        super(EmbedShell,self).__init__(
            config=config, ipython_dir=ipython_dir, user_ns=user_ns,
            user_module=user_module, custom_exceptions=custom_exceptions,
            usage=usage, banner1=banner1, banner2=banner2,
            display_banner=display_banner
        )

        self.exit_msg = exit_msg

        # don't use the ipython crash handler so that user exceptions aren't
        # trapped
        sys.excepthook = ultratb.FormattedTB(color_scheme=self.colors,
                                             mode=self.xmode,
                                             call_pdb=self.pdb)

    def init_sys_modules(self):
        pass

    def init_magics(self):
        super(EmbedShell, self).init_magics()
        self.register_magics(EmbeddedMagics)

    def __call__(self, header='', local_ns=None, module=None, dummy=None,
                 stack_depth=1, global_ns=None):
        """Activate the interactive interpreter.

        __call__(self,header='',local_ns=None,module=None,dummy=None) -> Start
        the interpreter shell with the given local and global namespaces, and
        optionally print a header string at startup.

        The shell can be globally activated/deactivated using the
        dummy_mode attribute. This allows you to turn off a shell used
        for debugging globally.

        However, *each* time you call the shell you can override the current
        state of dummy_mode with the optional keyword parameter 'dummy'. For
        example, if you set dummy mode on with IPShell.dummy_mode = True, you
        can still have a specific call work by making it as IPShell(dummy=False).
        """

        # If the user has turned it off, go away
        if not self.embedded_active:
            return

        # Normal exits from interactive mode set this flag, so the shell can't
        # re-enter (it checks this variable at the start of interactive mode).
        self.exit_now = False

        # Allow the dummy parameter to override the global __dummy_mode
        if dummy or (dummy != 0 and self.dummy_mode):
            return

        if self.has_readline:
            self.set_readline_completer()

        # self.banner is auto computed
        if header:
            self.old_banner2 = self.banner2
            self.banner2 = self.banner2 + '\n' + header + '\n'
        else:
            self.old_banner2 = ''

        # Call the embedding code with a stack depth of 1 so it can skip over
        # our call and get the original caller's namespaces.
        self.mainloop(local_ns, module, stack_depth=stack_depth, global_ns=global_ns)

        self.banner2 = self.old_banner2

        if self.exit_msg is not None:
            print self.exit_msg

    def mainloop(self, local_ns=None, module=None, stack_depth=0,
                 display_banner=None, global_ns=None):
        """Embeds IPython into a running python program.

        Input:

          - header: An optional header message can be specified.

          - local_ns, module: working local namespace (a dict) and module (a
          module or similar object). If given as None, they are automatically
          taken from the scope where the shell was called, so that
          program variables become visible.

          - stack_depth: specifies how many levels in the stack to go to
          looking for namespaces (when local_ns or module is None).  This
          allows an intermediate caller to make sure that this function gets
          the namespace from the intended level in the stack.  By default (0)
          it will get its locals and globals from the immediate caller.

        Warning: it's possible to use this in a program which is being run by
        IPython itself (via %run), but some funny things will happen (a few
        globals get overwritten). In the future this will be cleaned up, as
        there is no fundamental reason why it can't work perfectly."""
        
        if (global_ns is not None) and (module is None):
            class DummyMod(object):
                """A dummy module object for embedded IPython."""
                pass
            warnings.warn("global_ns is deprecated, use module instead.", DeprecationWarning)
            module = DummyMod()
            module.__dict__ = global_ns

        # Get locals and globals from caller
        if (local_ns is None or module is None) and self.default_user_namespaces:
            call_frame = sys._getframe(stack_depth).f_back

            if local_ns is None:
                local_ns = call_frame.f_locals
            if module is None:
                global_ns = call_frame.f_globals
                module = sys.modules[global_ns['__name__']]
        
        # Save original namespace and module so we can restore them after 
        # embedding; otherwise the shell doesn't shut down correctly.
        orig_user_module = self.user_module
        orig_user_ns = self.user_ns
        
        # Update namespaces and fire up interpreter
        
        # The global one is easy, we can just throw it in
        if module is not None:
            self.user_module = module

        # But the user/local one is tricky: ipython needs it to store internal
        # data, but we also need the locals. We'll throw our hidden variables
        # like _ih and get_ipython() into the local namespace, but delete them
        # later.
        if local_ns is not None:
            self.user_ns = local_ns
            self.init_user_ns()

        # Patch for global embedding to make sure that things don't overwrite
        # user globals accidentally. Thanks to Richard <rxe@renre-europe.com>
        # FIXME. Test this a bit more carefully (the if.. is new)
        # N.B. This can't now ever be called. Not sure what it was for.
        # And now, since it wasn't called in the previous version, I'm
        # commenting out these lines so they can't be called with my new changes
        # --TK, 2011-12-10
        #if local_ns is None and module is None:
        #    self.user_global_ns.update(__main__.__dict__)

        # make sure the tab-completer has the correct frame information, so it
        # actually completes using the frame's locals/globals
        self.set_completer_frame()

        with nested(self.builtin_trap, self.display_trap):
            self.interact(display_banner=display_banner)
        
        # now, purge out the local namespace of IPython's hidden variables.
        if local_ns is not None:
            for name in self.user_ns_hidden:
                local_ns.pop(name, None)
        
        # Restore original namespace so shell can shut down when we exit.
        self.user_module = orig_user_module
        self.user_ns = orig_user_ns

_embedded_shell = None

cfg = Config()
# new prompts
prompt_config = cfg.PromptManager
prompt_config.in_template = 'OPEN [\\#]: '
prompt_config.in2_template = '   .\\D.: '
prompt_config.out_template = 'Out<\\#>: '
