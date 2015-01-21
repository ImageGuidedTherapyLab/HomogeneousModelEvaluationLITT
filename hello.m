function hello
% This is the hello, world function written in M code

% Copyright 1997 The MathWorks, Inc.
% $Revision: 1.1.6.1 $
%
        %inputfilename  = sys.argl[1]
        %outputfilename = sys.argl[2]


        inputfilename  = 'optpp.mat.in.1'
        outputfilename = 'optpp.out.1'

        inputdata  = load(inputfilename )
        fprintf(1,'Hello, World\n' );
       
        %fncvalue = callsam(inputdata)
        fncvalue = 20.

        open (outputfilename )
        fprintf(1,fncvalue  );
