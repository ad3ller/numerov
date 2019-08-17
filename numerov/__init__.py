#! python
# -*- coding: utf-8 -*-

try:
    from .cy.core import radial_wf, radial_integral
except ImportError:
    from .core import radial_wf, radial_integral
except:
    raise
