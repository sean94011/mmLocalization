#include <AccelStepper.h>
#include <MultiStepper.h>

// TODO: Change 
// const long stepsPerRevolution=6400*8;  // only determines the torque and affect the speed
// // recommend to set it to 4 or 8 times the Microstep number 
// // const int NumPoint = 120;
// const int MaxNumPoint = 200;
// // const int D = 1248; // 1/2*lambda (1.95mm for 77GHz)
// const int D = 1600; // 1/2*lambda (2.5mm for 77GHz)
// const double spd = 10;

// Initialize the two steppers
AccelStepper xStepper(AccelStepper::DRIVER, 6, 7);
AccelStepper yStepper(AccelStepper::DRIVER, 8, 9);

// Create a MultiStepper Group
MultiStepper steppers;

//Record coordinate
long positions[2];
  
// Stepper xStepper(stepsPerRevolution, 6, 7);
// Stepper yStepper(stepsPerRevolution, 8, 9);

void setup() {
    //Currently assume mm/s speed, cm position
    //Multiplier will be different depending on motion stage used 
    //See steps/turn and then the drive length of on turn on linear actuator to define
    //Run to new position will drive motor until it reaches the location given (input is location, not distance to move)
    // xStepper.setSpeed(75);
    // yStepper.setSpeed(75);
    // delay(20);

    // Initialize the Serial Port
    Serial.begin(9600);

    // Configure each stepper
    xStepper.setMaxSpeed(100);
    yStepper.setMaxSpeed(100);

    // Add the two steppers to the MultiStpper Group to manage
    steppers.addStepper(xStepper);
    steppers.addStepper(yStepper);

    Serial.println("------------------------------------------");
    // Serial.println("-------------- x-plotter (max step 200)-----------------");
    // Serial.println("Enter Number of step (1/2*lambda step = 2.5mm)");
    // Serial.println("Positive: Move right,  Negative: Move left");
    Serial.println("Enter the 2D Position in [cm]");
    Serial.println("------------------------------------------");

}

void loop() {
  while (Serial.available() > 0) {
    // Parse the user input coordinate
    positions[0] = Serial.parseInt();
    positions[1] = Serial.parseInt();

    // Print the User Input
    printPosition(positions[0], positions[1]);

    // Set the destinations for the steppers
    steppers.moveTo(positions);

    // Blocks until all are in position
    steppers.runSpeedToPosition();

    // Delay to stablize
    delay(1000);
//     int NumStep[2] = {position2step(x) * (-1), position2step(y)};
//     maxNumStep = max(abs(NumStep[0]),abs(NumStep[1]));
//     xSign = signCheck(NumStep[0]);
//     ySign = signCheck(NumStep[1]);
    
//     if (Serial.read() == '\n') {
//       if (NumStep[0] > MaxNumPoint || NumStep[1] > MaxNumPoint){
//           Serial.println("Number of step exceeds max number of step.");
//           return;
//       }
//       else{
//         if(NumStep[0] >= 0 && NumStep[1] >= 0){
//             for (int k=0; k<NumStep; k++){
//               xStepper.step(D); // 1*lambda move left 
//             }                     
//         }
//         else{
//             for (int k=0; k<-NumStep; k++){
//               xStepper.step(-D); // 1*lambda move left 
//             }                               
//         }          
//       }
//     }    
  }
}

// total time = 4992*60/(20*6400) seconds

// Function: signCheck (inline)
// Description: check the moving direction of the steps
// Input: NumStep (int)
// Output: +1 or -1 (int)
inline int signCheck(int NumStep) {
  return (NumStep >= 0) ? 1 : -1;
}

// Function: smallerStep (inline)
// Description: check which direction has smaller number of steps
// Input: xNumStep (int), yNumStep (int)
// Output: 

inline int smallerStep(int xNumStep, int yNumStep)

// Function: position2step
// Description: Convert coordinate to number of steps
// Input: Single coordinate value (float)
// Ouput: Number of Steps (int)
int position2step(float coord) {
    // TODO: Get the mathematical formula of the conversion
    int NumStep = coord; // temporary return value
    return NumStep;
}

// TODO
// Function: reset
// Description: Reset the motion stage to the motor's original state
// Input: ?
// Output: N/A
void reset() {
    
}

// Function: Print Position (Inline Function)
// Description: Print the Input Coordinate
// Input: x coordinate (int), y coordinate (int)
// Output: N/A
inline void printPosition(int x, int y) {
    char msg[20];
    sprintf(msg, "(x, y) = (%d,%d)", x, y);
    Serial.println(msg);
}