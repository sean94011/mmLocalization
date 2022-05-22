%% Initialize Serialport Connection

%xStepper = arduino('/dev/tty.usbmodem21401','Uno');
%yStepper = arduino('/dev/tty.usbmodem21301','Uno');

delete(instrfindall);
SPEED_CONTROL=100;
DIRECTION=-1;
%DIRECTION=1;
xStepper = serialport("/dev/cu.usbmodem21401",9600)
yStepper = serialport("/dev/cu.usbmodem21301",9600)
configureTerminator(xStepper,"CR");
configureTerminator(yStepper,"CR");
flush(xStepper);
flush(yStepper);
