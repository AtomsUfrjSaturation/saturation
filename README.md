### This project simplifies the process of creating a Hysys extension with custom logic. 

The repo's root contains C# code responsible for communicating with Hysys using Windows' COM 
(you need the Aspen Developer Kit for the interop DLL `Interop.HYSYS.1.7.dll`). 

Inside `kernel` there is a PC-SAFT modelling logic written in Rust that 
leverages FFI for exchanging the values (such as temperature, pressure etc) between the two 
languages. 

#### Note: receiving a stream's components list
This data is defined at runtime, so it's available only
when the user is interacting with the extension in Hysys.
To access it: 

```csharp
try {
    dynamic components = Feed.Flowsheet.FluidPackage.Components;
    string componentNames = string.Join(", ", components.Names);
    // remaining code...
} catch (Exception e) {
    Logger(e.ToString());
}
```