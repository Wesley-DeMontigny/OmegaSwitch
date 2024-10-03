#ifndef MOVE_HPP
#define MOVE_HPP

class ParameterNode;

class Move{
    public:
        virtual double update()=0;
        virtual void tune()=0;
        void markAccepted() {acceptedSinceTune++; acceptedCount++; countSinceTune++;}
        void markRejected() {countSinceTune++; rejectedCount++;}
        void clearRecord() {countSinceTune = 0; acceptedCount = 0; acceptedSinceTune = 0; rejectedCount = 0;}
    protected:
        int countSinceTune = 0;
        int acceptedSinceTune = 0;
        int rejectedCount = 0;
        int acceptedCount = 0;
};

#endif