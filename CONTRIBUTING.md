Contributions are welcomed via merge requests on Gitlab.com.
First contact the developers prior to beginning and propose your ideas for a new feature or bug fix.
Then implement your code.

Submit a merge request on gitlab. Multiple developers and/or users will review requested changes and make comments.
The lead developer will merge into the mainline after the review is complete and approved.

The rest of this file will be used as a checklist to review the merge request.

# Features

## Implement functionality in a general and flexible fashion

SOMA intends to provide flexibility for the user. Your merge request should provide something that is applicable
to a variety of use-cases and not just the one thing you might need it to do for your research. Speak to the lead
developers before writing your code, and they will help you make design choices that allow flexibility.

## Do not degrade performance of existing code paths

New functionalities should only activate expensive code paths when they are requested by the user. Do not slow down
existing code.

## Optimize for the current accelerator generation

Use OpenACC/OpenMP if your feature is compute intensive. Test optimize your code for CPU execution and accelerators/GPUs.

# Version control

## Base your work off the correct branch

Bug fixes should be based on `maint`. New features should be based on `master`.

## Propose a single set of related changes

Changes proposed in a single topic branch / merge request should all be related to each other. Don't propose too
many changes at once, review becomes challenging. Multiple new features that are loosely coupled should be completed
in separate topic branches. It is OK if the branch for `feature2` is based on `feature1` - as long as it is made clear
that `feature1` should be merged before the review of `feature2`. It is better to merge both `feature1` and `feature2`
into a temporary integration branch during testing.

## Keep changes to a minimum

Don't go and "fix" spelling errors all over the code, or make lots of whitespace changes along with a new feature.
If there are spelling errors to fix, propose that in a separate merge request.

## Agree to the contributor agreement

All contributors must agree to the Contributor Agreement ([ContributorAgreement](ContributorAgreement)) before their merge request can be merged.

# Source code

## Use our style

It is important to have a consistent style throughout the source code.
For this matter we use a fixed style. The exact style can be forced with a script reformating the source code: c_src/force_indention.sh.
The GNU indent tool is required to reformat the source code. Run the script before you push your contribution. If you added new source code and they are not automaticaly detected, add the files manually to the script.

Readibilty of your code matters.

## Document code with comments

Use doxygen header comments for structs, functions, etc... Also comment complex sections of code so that other
developers can understand them.

## Compiles without warnings

Your changes should compile without warnings.

# Tests

## Write  tests

All new functionality in SOMA should be tested with automatic tests.
Existing test should pass with the merge request.

## Validity tests

In addition to the automatic tests, the developer should run research-scale simulations using the new functionality and
ensure that it behaves as intended.

## Add developer to the authors list

Developers need to be credited for their work. In addition to the authors list add your name to the changes in the  change log. 
contributed to the code.

## Propose a change log

Propose a short concise entry describing the change for use in the change log.

## Reference

Based on the CONTRIBUTING.md of the hoomd-blue project: https://bitbucket.org/glotzer/hoomd-blue
