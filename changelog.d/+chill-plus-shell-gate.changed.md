Users who run chill-plus on a cutoff list now see different counts.
This is a bug fix in the shell count.
The old code counted bonds after the list had been cut to four.
A crystal row needs four neighbours inside the cutoff.
On the mixed TIP4P fixture, interClathrate moves from 12 to 9.
Water on that fixture moves from 238 to 241.
Default q3 uses up to four nearest neighbours in the cutoff.
Interfacial rows keep the neighbour test from the published table.
Users who want the OVITO labels use the ovito-chill+ command.
That command averages q3 on the full shell and skips the test.
Plain CHILL keeps four closest neighbours and Moore's test.
A reader error now returns exit code 2 instead of aborting the process.
We print reCubic and reHex as plain text.
A path that ends in xyz and cannot be opened returns an error.
An exception outside std::exception returns exit code 2.
