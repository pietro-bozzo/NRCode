function ids = printRegs_(session)

R = regions(session,verbose=false);
ids = R.ids;