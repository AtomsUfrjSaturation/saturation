### A fast and modern saturation calculation plugin for Aspen Hysys

#### Overview
This project delivers a saturation calculation plugin for Hysys. It offers a modern build approach to accurately determine water saturation levels within the Hysys environment. Additionally, the structure and design of this plugin serve as a secondary benefit, simplifying the creation of Hysys extensions with custom logic.

#### Features
- **Saturation calculations**: At its core, this plugin is made for fast and reliable saturation calculations, particularly for water.

- **Custom logic integration**: While the primary focus is on saturation calculations, the plugin also demonstrates how custom logic can be seamlessly integrated into Hysys extensions, offering a valuable resource for developers.

#### Repository structure
- **C# Interface with Hysys** (Root directory): Contains C# code for interfacing with Hysys using Windows' COM. This allows for effective communication between the Hysys software and the plugin's custom logic.

- **PC-SAFT Modelling in Rust** (`kernel` directory): The kernel directory houses the PC-SAFT modelling logic, implemented in Rust. It employs Foreign Function Interface (FFI) for data exchange (like temperature, pressure, etc.) between C# and Rust.

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
