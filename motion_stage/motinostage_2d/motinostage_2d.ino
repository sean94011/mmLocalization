#include <AccelStepper.h>
#include <MultiStepper.h>

// Function: Print Position (Inline Function)
// Description: Print the Input Coordinate
// Input: x coordinate (int), y coordinate (int)
// Output: N/A
inline void printPosition(int x, int y) {
    char msg[20];
    sprintf(msg, "(x, y) = (%d,%d)", x, y);
    Serial.println(msg);
}

// Initialize the two steppers
const int xPulPin = 6;
const int xDirPin = 7;
const int yPulPin = 8;
const int yDirPin = 9;
AccelStepper xStepper(AccelStepper::DRIVER, xPulPin, xDirPin);
AccelStepper yStepper(AccelStepper::DRIVER, yPulPin, yDirPin);

// Create a MultiStepper Group
MultiStepper steppers;

// Record coordinate
long positions[2] = {0,0};

// Record Speed
float speed = 0;

void setup() {
    // Initialize the Serial Port
    Serial.begin(9600);

    // Add the two steppers to the MultiStpper Group to manage
    steppers.addStepper(xStepper);
    steppers.addStepper(yStepper);

    Serial.println("------------------------------------------");
    Serial.println("Enter 'speed, x, y' [mm/s,mm,mm]");
    Serial.println("------------------------------------------");
}

void loop() {
  while (Serial.available() > 0) {
    // Parse the user input coordinate
    speed = Serial.parseFloat();
    positions[0] = round(Serial.parseFloat());
    positions[1] = round(Serial.parseFloat());

    // Print the User Input
    printPosition(positions[0], positions[1]);

    // Configure each stepper
    xStepper.setMaxSpeed(round(speed));
    yStepper.setMaxSpeed(round(speed));

    // Set the destinations for the steppers
    positions[0] *= -1;
    steppers.moveTo(positions);

    // Blocks until all are in position
    steppers.runSpeedToPosition();

    // Delay to stablize
    delay(1000);
  }
}

