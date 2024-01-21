### Simplifying the creation of Hysys extensions with custom logic 

This project focuses on streamlining the process of developing Hysys extensions that incorporate custom logic. 

The repository's root directory contains C# code, which is responsible for interfacing with Hysys via Windows' COM. This setup facilitates communication between the Hysys software and our custom logic.

Within the `kernel` directory, you'll find the PC-SAFT modelling logic, implemented in Rust. This component uses Foreign Function Interface (FFI) to facilitate data exchange (such as temperature, pressure, etc.) between C# and Rust.

#### Note 0: Handling File Mode Changes
Working across different operating systems (Unix and Windows) often leads to file mode changes. These changes can obscure actual code modifications when comparing environments. To mitigate this issue, execute the following command at the root of the repository to ignore such changes:
```
git config core.fileMode false
``` 

#### Note 1: Accessing a Stream's Components List
This data is dynamically defined at runtime and is only accessible when the user interacts with the extension in Hysys. To retrieve the list of components in a stream, use the following C# code snippet:

```csharp
try {
    dynamic components = Feed.Flowsheet.FluidPackage.Components;
    string componentNames = string.Join(", ", components.Names);
    // Add additional code here as needed...
} catch (Exception e) {
    Logger(e.ToString());
}
```
