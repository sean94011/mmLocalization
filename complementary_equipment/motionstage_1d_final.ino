#include <AccelStepper.h>
#include <MultiStepper.h>

// Function: Print Position (Inline Function)
// Description: Print the Input Coordinate
// Input: x coordinate (int), y coordinate (int)
// Output: N/A
inline void printPosition(int t, int x, int y, int y2) {
    char msg[20];
    sprintf(msg, "time = %d, (x, y, y2) = (%d,%d,%d)", t, x, y, y2);
    Serial.println(msg);
}

// Initialize the two steppers
const int xPulPin = 6;
const int xDirPin = 7;
const int yPulPin = 8;
const int yDirPin = 9;
// Initialize the second motionstage
const int y2PulPin = 4;
const int y2DirPin = 5;

AccelStepper xStepper(AccelStepper::DRIVER, xPulPin, xDirPin);
AccelStepper yStepper(AccelStepper::DRIVER, yPulPin, yDirPin);
AccelStepper y2Stepper(AccelStepper::DRIVER, y2PulPin, y2DirPin);

// Create a MultiStepper Group
MultiStepper steppers;

// Record coordinate
long pos[3] = {0,0,0};

// Record Speed & Time
float speed = 0;
float time = 0;

// Record absolute position of motionstage relative to it's origin
float curx = 0; // Movement in the left and right direction (towards the first radar)
float cury = 0; // Movement in the forward and backward direction (towards the second radar)
float cury2 = 0; // Movement in the forward and backward direction (towards the second radar)

// Function: Take user input for coordinates
// Description: Receive coordinates from serial input
// Input: running time (int), x coordinate (int), y coordinate (int), y2 coordinate for second motionstage (int)
// Output: N/A
void setup() {
    // Initialize the Serial Port
    Serial.begin(9600);

    // Add the two steppers to the MultiStpper Group to manage
    steppers.addStepper(xStepper);
    steppers.addStepper(yStepper);
    steppers.addStepper(y2Stepper);

    Serial.println("------------------------------------------");
    Serial.println("Enter 'time, x, y, y2' [sec,mm,mm,mm]");
    Serial.println("------------------------------------------");
}

// Function: Move the motionstage
// Description: Move the motionstage to the required x,y,y2 coordinates in the required amount of time
// Input: N/A
// Output: N/A (Motion)
void loop() {
  while (Serial.available() > 0) {
    // Parse the user input coordinate
    time = Serial.parseFloat();
    
    pos[0] = round(Serial.parseFloat()*160);
    pos[1] = round(Serial.parseFloat()*160);
    pos[2] = round(Serial.parseFloat()*160);
    
    // Check the maximum distance
    int maxDist = max(max(abs(pos[0] - curx),abs(pos[1] - cury)),abs(pos[2] - cury2));

    // Compute the velocity by the designated time duration
    speed = maxDist / time;
    curx = pos[0];
    cury = pos[1];
    cury2 = pos[2];
    
    // Print the User Input
    printPosition(time, pos[0], pos[1], pos[2]);

    // Configure each stepper
    xStepper.setMaxSpeed(round(speed));
    yStepper.setMaxSpeed(round(speed));
    y2Stepper.setMaxSpeed(round(speed));
    
    // Set the destinations for the steppers
    pos[0] *= -1; //change to pos if want to reverse xy directions
    steppers.moveTo(pos);
    
    // Blocks until all are in position
    steppers.runSpeedToPosition();

    while (Serial.available() > 0) Serial.read();

    // Delay to stablize
    delay(1000);
  }
}
