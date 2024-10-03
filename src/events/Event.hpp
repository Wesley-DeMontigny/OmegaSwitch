#ifndef EVENT_HPP
#define EVENT_HPP

class Event{
    public:
        virtual void initialize()=0;
        virtual void call(int iteration)=0;
    protected:
        int interval;
};

#endif