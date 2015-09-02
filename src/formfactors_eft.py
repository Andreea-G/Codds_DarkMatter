"""
Copyright (c) 2015 Andreea Georgescu

Created on Sat Nov 22 17:20:30 2014

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import numpy as np
from globalfnc import *

""" Effective Field Theory (EFT) Form Factors FFSigmaPPJ and FFSigmaPJ
(See Appendix A.3. of Fitzpatrick et al http://arxiv.org/abs/1203.3542)
"""
FFSigmaPPJ = {(19, 9): np.array([[lambda y: 0.90278 + y * (-2.37144 + y * (2.3531 + y * (-1.04517 + 0.175359 * y))),
                                  lambda y: -0.0165506 + y * (0.050948 + y * (-0.0510308 + y * (0.0199287 - 0.00236734 * y)))],
                                 [lambda y: -0.0165506 + y * (0.050948 + y * (-0.0510308 + y * (0.0199287 - 0.00236734 * y))),
                                  lambda y: 0.000303421 + y * (-0.00107102 + y * (0.00114207 + y * (-0.000347594 + 0.000031959 * y)))]]),
              (23, 11): np.array([[lambda y: 0.136361 + y * (-0.266817 + y * (0.457598 + y * (-0.112141 + 0.00828489 * y))),
                                   lambda y: 0.0109688 + y * (-0.029996 + y * (0.0217192 + y * (-0.00897498 + 0.000591796 * y)))],
                                  [lambda y: 0.0109688 + y * (-0.029996 + y * (0.0217192 + y * (-0.00897498 + 0.000591796 * y))),
                                   lambda y: 0.000882318 + y * (-0.00309927 + y * (0.00398861 + y * (-0.00203246 + 0.00040932 * y)))]]),
              (73, 32): np.array([[lambda y: 0.000102041 + y * (-0.00096119 + y * (0.00335189 + y * (-0.00364089 + y * (0.0019991 + y * (-0.000462685 + y * (0.000040244 + y * (-3.45152e-7 + 1.2863e-9 * y))))))),
                                   lambda y: 0.00612884 + y * (-0.0386733 + y * (0.0575045 + y * (-0.0535467 + y * (0.0257882 + y * (-0.0065118 + y * (0.000816404 + y * (-0.0000402374 + 7.36917e-7 * y)))))))],
                                  [lambda y: 0.00612884 + y * (-0.0386733 + y * (0.0575045 + y * (-0.0535467 + y * (0.0257882 + y * (-0.0065118 + y * (0.000816404 + y * (-0.0000402374 + 7.36917e-7 * y))))))),
                                   lambda y: 0.368112 + y * (-1.17814 + y * (2.27237 + y * (-1.98163 + y * (1.01909 + y * (-0.296483 + y * (0.056608 + y * (-0.00599025 + 0.000951084 * y)))))))]]),
              (127, 53): np.array([[lambda y: 0.130131 + y * (-0.488539 + y * (1.83369 + y * (-2.79561 + y * (2.70421 + y * (-1.58275 + y * (0.528399 + y * (-0.0920503 + y * (0.00666861 + y * (-6.28403e-12 + 1.3501e-20 * y))))))))),
                                    lambda y: 0.0322868 + y * (-0.125179 + y * (0.26445 + y * (-0.296832 + y * (0.214155 + y * (-0.0982407 + y * (0.0271252 + y * (-0.00422686 + y * (0.000330271 + y * (-9.72551e-6 + 1.02434e-14 * y)))))))))],
                                   [lambda y: 0.0322868 + y * (-0.125179 + y * (0.26445 + y * (-0.296832 + y * (0.214155 + y * (-0.0982407 + y * (0.0271252 + y * (-0.00422686 + y * (0.000330271 + y * (-9.72551e-6 + 1.02434e-14 * y))))))))),
                                    lambda y: 0.00801063 + y * (-0.0320427 + y * (0.0534333 + y * (-0.0461654 + y * (0.0247382 + y * (-0.00858068 + y * (0.00190294 + y * (-0.000262031 + y * (0.0000213702 + y * (-9.3627e-7 + 1.68751e-8 * y)))))))))]]),
              (129, 54): np.array([[lambda y: 0.000210575 + y * (-0.00150295 + y * (0.00431238 + y * (-0.00630033 + y * (0.00491101 + y * (-0.00199528 + y * (0.00041825 + y * (-0.0000422967 + y * (1.6268e-6 + y * (3.31573e-14 + y * 1.68953e-22))))))))),
                                    lambda y: 0.00720992 + y * (-0.0380102 + y * (0.0824299 + y * (-0.0971526 + y * (0.067932 + y * (-0.0271548 + y * (0.00613825 + y * (-0.000734195 + y * (0.0000354172 + y * (-7.27679e-8 + y * -7.41579e-16)))))))))],
                                   [lambda y: 0.00720992 + y * (-0.0380102 + y * (0.0824299 + y * (-0.0971526 + y * (0.067932 + y * (-0.0271548 + y * (0.00613825 + y * (-0.000734195 + y * (0.0000354172 + y * (-7.27679e-8 + y * -7.41579e-16))))))))),
                                    lambda y: 0.246862 + y * (-0.840937 + y * (1.4482 + y * (-1.46722 + y * (0.944902 + y * (-0.37255 + y * (0.0890997 + y * (-0.0120716 + y * (0.000755735 + y * (-3.08386e-6 + y * 3.255e-9)))))))))]]),
              (131, 54): np.array([[lambda y: 0.0000578995 + y * (-0.000032744 + y * (0.000127388 + y * (-0.000624606 + y * (0.000875357 + y * (-0.000528633 + y * (0.000151948 + y * (-0.0000197785 + y * (9.38122e-7 + y * (-1.04681e-13 + y * 5.94262e-21))))))))),
                                    lambda y: 0.00225524 + y * (0.00315534 + y * (-0.0106065 + y * (-0.000772228 + y * (0.0192839 + y * (-0.0176928 + y * (0.00663425 + y * (-0.00112357 + y * (0.0000683327 + y * (-3.76191e-8 + y * 3.59162e-15)))))))))],
                                   [lambda y: 0.00225524 + y * (0.00315534 + y * (-0.0106065 + y * (-0.000772228 + y * (0.0192839 + y * (-0.0176928 + y * (0.00663425 + y * (-0.00112357 + y * (0.0000683327 + y * (-3.76191e-8 + y * 3.59162e-15))))))))),
                                    lambda y: 0.0878436 + y * (0.295485 + y * (-0.232083 + y * (-0.4651 + y * (1.24324 + y * (-1.07229 + y * (0.438103 + y * (-0.085951 + y * (0.006652 + y * (-7.67519e-6 + y * 2.24565e-9)))))))))]]),
              }

FFSigmaPJ = {(19, 9): np.array([[lambda y: 1.80556 + y * (-4.8508 + y * (4.87922 + y * (-2.17774 + 0.363912 * y))),
                                 lambda y: -0.0331012 + y * (0.0814569 + y * (-0.0511344 + y * (-0.00141634 + 0.00602374 * y)))],
                                [lambda y: -0.0331012 + y * (0.0814569 + y * (-0.0511344 + y * (-0.00141634 + 0.00602374 * y))),
                                 lambda y: 0.000606843 + y * (-0.00135635 + y * (0.000265925 + y * (0.000549795 + 0.0000997092 * y)))]]),
             (23, 11): np.array([[lambda y: 0.272723 + y * (-0.824074 + y * (1.18533 + y * (-0.476888 + 0.059333 * y))),
                                  lambda y: 0.0219376 + y * (-0.0577543 + y * (0.0359879 + y * (-0.00300412 - 0.000363159 * y)))],
                                 [lambda y: 0.0219376 + y * (-0.0577543 + y * (0.0359879 + y * (-0.00300412 - 0.000363159 * y))),
                                  lambda y: 0.00176464 + y * (-0.00395928 + y * (0.00228214 + y * (0.0000194667 + 2.36537e-6 * y)))]]),
             (73, 32): np.array([[lambda y: 0.000204083 + y * (-0.000463558 + y * (0.00148348 + y * (-0.00229056 + y * (0.00205049 + y * (-0.000756444 + y * (0.000102527 + y * (-2.05037e-6 + 3.00769e-8 * y))))))),
                                  lambda y: 0.0122577 + y * (-0.0531362 + y * (0.0723673 + y * (-0.0669153 + y * (0.038403 + y * (-0.0119103 + y * (0.00184147 + y * (-0.000139208 + 4.97489e-6 * y)))))))],
                                 [lambda y: 0.0122577 + y * (-0.0531362 + y * (0.0723673 + y * (-0.0669153 + y * (0.038403 + y * (-0.0119103 + y * (0.00184147 + y * (-0.000139208 + 4.97489e-6 * y))))))),
                                  lambda y: 0.736223 + y * (-4.71068 + y * (11.6726 + y * (-12.5199 + y * (6.91299 + y * (-2.0712 + y * (0.347982 + y * (-0.030887 + 0.00191149 * y)))))))]]),
             (127, 53): np.array([[lambda y: 0.260263 + y * (-1.59356 + y * (5.29612 + y * (-8.89098 + y * (8.70123 + y * (-4.9109 + y * (1.53875 + y * (-0.245616 + y * (0.0157491 + y * (-1.18163e-10 + 9.23323e-19 * y))))))))),
                                   lambda y: 0.0645735 + y * (-0.455982 + y * (1.28087 + y * (-1.79425 + y * (1.41544 + y * (-0.651004 + y * (0.174051 + y * (-0.0261935 + y * (0.00201904 + y * (-0.0000604279 + 5.34917e-13 * y)))))))))],
                                  [lambda y: 0.0645735 + y * (-0.455982 + y * (1.28087 + y * (-1.79425 + y * (1.41544 + y * (-0.651004 + y * (0.174051 + y * (-0.0261935 + y * (0.00201904 + y * (-0.0000604279 + 5.34917e-13 * y))))))))),
                                   lambda y: 0.0160213 + y * (-0.12817 + y * (0.370306 + y * (-0.482942 + y * (0.335727 + y * (-0.135075 + y * (0.0327184 + y * (-0.00480083 + y * (0.000414369 + y * (-0.0000192222 + 3.67902e-7 * y)))))))))]]),
             (129, 54): np.array([[lambda y: 0.00042115 + y * (-0.00186625 + y * (0.00648086 + y * (-0.0119809 + y * (0.0166818 + y * (-0.012071 + y * (0.00413563 + y * (-0.000626931 + y * (0.0000341264 + y * (1.28862e-12 + 1.21646e-20 * y))))))))),
                                   lambda y: 0.0144198 + y * (-0.0796792 + y * (0.230627 + y * (-0.414422 + y * (0.427739 + y * (-0.237843 + y * (0.0697309 + y * (-0.0101269 + y * (0.000584933 + y * (-2.82804e-6 - 5.33937e-14 * y)))))))))],
                                  [lambda y: 0.0144198 + y * (-0.0796792 + y * (0.230627 + y * (-0.414422 + y * (0.427739 + y * (-0.237843 + y * (0.0697309 + y * (-0.0101269 + y * (0.000584933 + y * (-2.82804e-6 - 5.33937e-14 * y))))))))),
                                   lambda y: 0.493723 + y * (-3.26846 + y * (8.78659 + y * (-12.3673 + y * (9.84395 + y * (-4.51071 + y * (1.1774 + y * (-0.16417 + y * (0.00997399 + y * (-0.0000926412 + 2.3436e-7 * y)))))))))]]),
             (131, 54): np.array([[lambda y: 0.000115799 + y * (-0.000893646 + y * (0.00154676 + y * (0.00149625 + y * (-0.000692549 + y * (-0.00118632 + y * (0.00080076 + y * (-0.0001592 + y * (0.0000101291 + y * (-3.38696e-12 + 2.9816e-19 * y))))))))),
                                   lambda y: 0.00451048 + y * (-0.0385549 + y * (0.0947553 + y * (-0.0375968 + y * (-0.0774151 + y * (0.0873457 + y * (-0.0346308 + y * (0.00586233 + y * (-0.000352003 + y * (-9.39713e-7 + 1.6094e-13 * y)))))))))],
                                  [lambda y: 0.00451048 + y * (-0.0385549 + y * (0.0947553 + y * (-0.0375968 + y * (-0.0774151 + y * (0.0873457 + y * (-0.0346308 + y * (0.00586233 + y * (-0.000352003 + y * (-9.39713e-7 + 1.6094e-13 * y))))))))),
                                   lambda y: 0.175687 + y * (-1.64768 + y * (5.75397 + y * (-9.70477 + y * (9.11141 + y * (-4.85381 + y * (1.43134 + y * (-0.214356 + y * (0.0124188 + y * (0.0000657353 + 8.81625e-8 * y)))))))))]]),
             }


def FFElementQ(Z):
    """ Checks if the EFT Form Factors above have been implemented for the element
    with atomic number Z.
    Input:
        Z: ndarray
            List or atomic numbers.
    Returns:
        1 if yes and 0 if not.
    """
    check = [math.trunc(z) in np.array([9, 11, 32, 53, 54]) for z in Z]
    return np.array([1 if c else 0 for c in check])
