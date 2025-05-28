#pragma once 
#include <iostream>

template<typename T>
class ArbitrageError: public std::exception 
{
    public:
        const char* what() const noexcept override 
        {    
            if (_cachedMessage.empty()) {
                _cachedMessage = get_message();  // Safe: store in member
            }
            return _cachedMessage.c_str();
        };

        T get_code() const {return static_cast<T>(code_);};
        explicit ArbitrageError(const T code): code_(code){};
        virtual ~ArbitrageError() = default;
    protected: 
        
        std::string get_message() const 
        {
            return "Error [" + get_branch() + "] [" +  get_module_name() + "] [" + std::to_string(static_cast<int>(code_)) + "] : " + get_error_message();
        };
        virtual std::string get_branch() const = 0;
        virtual std::string get_module_name() const = 0; 
        virtual std::string get_error_message() const = 0; 
    private:
        const T code_;
        mutable std::string _cachedMessage;  // must be mutable to modify in const what()
};

template<typename T>
class CoreArbitrageError: public ArbitrageError<T>
{
    public:
        explicit CoreArbitrageError(const T code): ArbitrageError<T>(code){};
        virtual ~CoreArbitrageError() = default;
    protected: 
        std::string get_branch() const override {return "core";}
        virtual std::string get_module_name() const override = 0; 
        virtual std::string get_error_message() const override = 0; 
};



