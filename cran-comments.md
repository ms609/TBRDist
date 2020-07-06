This release should address the problem identified by clang-ASAN in CRAN's 
post-acceptance tests (see e-mail from Prof Ripley, 2020-07-01).

I am not aware of (and r-packages-devel has not suggested) any way for me to 
test this myself, as a Windows user, as the error is only identified by
Fedora-clang-asan, which is not available via rhub.


r-oldrel-osx-x86_64 also shows an installed size > 5Mb, due to a libs 
sub-directory of 5.5 Mb; this likely reflects the use of Boost headers,
and I'm not sure whether anything can be done about this.


## Test environments
* Windows 10 on local machine, R 4.0.2
* ubuntu 16.04.6 LTS (on travis-ci), R 3.4.0, release and devel, via [Travis CI](https://travis-ci.org/ms609/TBRDist)
* Mac OS X 10.13.6 (on travis-ci), R release
* win-builder, with `check_win_devel()`, R devel
* R-hub, with `check_for_cran()` and `check_with_sanitizers()`

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies for this package.
