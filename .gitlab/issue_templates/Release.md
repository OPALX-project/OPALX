**source code and binary**
  * [ ] create branch YEAR.N
  * [ ] create issue "Release version YEAR.N" and merge request
  * [ ] update version string in Doxyfile
  * [ ] update the version string in CMakeLists.txt and commit
  * [ ] wait for approval of MR and merge
  * [ ] tag version YEAR.N.0
  * [ ] upload source tar-ball to `/afs/psi.ch/project/amas/webhosting/Downloads/OPAL/src`
  * [ ] compile new binary for Linux
  * [ ] upload Linux binary package to `/afs/psi.ch/project/amas/webhosting/Downloads/OPAL/package/`
  * [ ] compile new binary for macOS
  * [ ] upload macOS binary package to `/afs/psi.ch/project/amas/webhosting/Downloads/OPAL/package/`

**manual/documentation**
  * [ ] setup a new branch for the new version of the manual
  * [ ] fix version, branches and links in `Manual.attributes`
  * [ ] clone repository into `/afs/psi.ch/project/amas/webhosting/opal/Documentation/x.y` and checkout new branch (`git clone https://gitlab.psi.ch/OPAL/documentation/manual.git`)
  * [ ] add links to the binaries in the wiki
  * [ ] update https://gitlab.psi.ch/OPAL/src/wikis/For-Developers/Compile-OPAL
  * [ ] compile the change log/release notes and publish it in the wiki: https://gitlab.psi.ch/OPAL/src/wikis/ReleaseNotes
  * [ ] review the file `src/addToDoxygenMainPage.h`
  * [ ] build Doxygen documentation
  * [ ] update https://gitlab.psi.ch/OPAL/src/wikis/home
  * [ ] update https://gitlab.psi.ch/OPAL/src/wikis/regression-tests

**tracker**
  * [ ] new milestone for `OPAL x.(y+1)`
  * [ ] update labels and milestones in issues

**regression-tests**

  * [ ] create new branch x.y
  * [ ] setup the regression-tests to run the new version on opalrunner.psi.ch

**varia**
  * [ ] PSI module
  * [ ] write e-mail to mailing list
